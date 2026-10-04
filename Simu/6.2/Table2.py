# ==============================================================================
# Script: Table2.py
# Description: Simulation study for the Multiple Predictor setting (Section 6.2).
#              Fits MTDR to two-dimensional distributions with one reference
#              and two random predictors, using input-convex neural networks
#              (ICNNs) and a fixed-support entropic barycenter approximation.
#
# Method: Multi-transport Distributional Regression (MTDR)
# Key Algorithm: Block-coordinate optimization of three maps and simplex weights.
#
# Output:
#   1. --out_csv: Trial-level results, updated after each successful replication.
#   2. --summary_csv: Monte Carlo means and sample SDs, appended across runs.
#      test_sinkhorn_loss_mean/std support the prediction comparison in Table 2.
#      The script does not assemble or round the manuscript table.
#
# Data:
#   Synthetic data are generated within each replication; no data files needed.
#   The reference is a truncated N(0, I2) law. The two random predictors are
#   independent draws from the same distribution-generating mechanism; they
#   are not identical realized measures. Observed and hidden particle clouds
#   are sampled independently by default. Only training responses are noisy.
#
# Dependencies:
#   Python 3, numpy, pandas, scipy, torch (PyTorch), geomloss.
#   CUDA is used when available; otherwise the script uses the CPU.
#
# Usage (run from this folder; use --help to list all options):
#   Quick execution check, NOT a full Monte Carlo experiment:
#     python Table2.py --N 10 --N_val 3 --N_test 3 --M 50 --trials 1 \
#         --outer_epochs 2 --inner_epochs 1 --alpha_epochs 1
#
#   Example Monte Carlo configuration (one weight/sample-size setting):
#     python Table2.py --N 50 --N_val 15 --N_test 15 --M 200 --trials 50 \
#         --alpha_true 0.3,0.35,0.35 --gpu_id 0 \
#         --out_csv p2_N50_M200.csv --summary_csv p2_summary.csv
#
# Parameter Conventions:
#   --alpha_true gives alpha0*, alpha1*, alpha2* in that order; nonnegative
#   inputs with positive sum are normalized. alpha0* is the reference weight.
#   --N is the training sample size n; --M is the particles per input measure m.
#   --N_val and --N_test are separate counts, not fractions of --N.
#   Replication t uses seed + 1000*t. The default --trials is 1; specify the
#   desired number of replications explicitly for a Monte Carlo summary.
#
# Evaluation and Numerical Notes:
#   - Barycenter computation uses squared Euclidean cost, default epsilon=0.25,
#     and 25 iterations on a 12-by-12 support grid over [-2.4,2.4]^2. This
#     numerical support differs from the input sampling domain [-2,2]^2.
#   - Training and evaluation use GeomLoss Sinkhorn divergence (p=2, blur=0.10).
#     test_sinkhorn_loss is its mean over clean test responses, without taking
#     a square root. It is not the exact matching-based W2^2 used in Table1.py.
#     epsilon_bary controls barycenter computation; blur_loss controls the
#     prediction discrepancy. They are distinct regularization settings.
#   - T0_L2/T1_L2/T2_L2 are grid-based root-mean-squared map discrepancies.
#     These and weight errors are diagnostics, not identifiability guarantees.
#   - Gaussian support noise and learned maps are not clipped. The anchor
#     penalty is soft; it does not impose a hard constraint on the maps.
#   - Early stopping and checkpoint selection use validation data, not test data.
#   - Only successful trials enter summaries; inspect trials_successful and logs.
#     A single successful trial has no sample SD (reported as a missing value).
#   - Output paths are relative to the working directory; parent folders must
#     exist. Use distinct raw CSV names for different settings/runs, since the
#     raw file is overwritten. Summary rows accumulate across runs.
#   - Seeds do not guarantee bitwise reproducibility across devices/versions.
#     This script is intended for CLI use: argument parsing and device setup
#     occur at module level.
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 1: Runtime Configuration and Command-Line Options
# ------------------------------------------------------------------------------

import os
import sys

# Avoid local code.py shadowing Python's standard-library code module during torch import.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path = [
    p for p in sys.path
    if os.path.abspath(p or os.getcwd()) != SCRIPT_DIR
]

# Limit low-level numerical-library threads before importing torch/numpy-heavy packages.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import copy
import math
import time

parser = argparse.ArgumentParser(
    description="p=2 MTDR with ICNN maps and fixed-support entropic Wasserstein barycenters"
)
parser.add_argument(
    "--N",
    type=int,
    default=50,
    help="Number of training distributional observations",
)
parser.add_argument(
    "--N_val",
    type=int,
    default=15,
    help="Number of validation distributional observations",
)
parser.add_argument(
    "--N_test",
    type=int,
    default=15,
    help="Number of test distributional observations",
)
parser.add_argument(
    "--M",
    type=int,
    default=200,
    help="Particles per predictor/reference distribution",
)
parser.add_argument("--trials", type=int, default=1, help="Number of Monte Carlo replications")
parser.add_argument(
    "--grid_n",
    type=int,
    default=12,
    help="Barycenter support grid side length; support size is grid_n^2",
)
parser.add_argument(
    "--support_low",
    type=float,
    default=-2.4,
    help="Lower bound of fixed barycenter support grid",
)
parser.add_argument(
    "--support_high",
    type=float,
    default=2.4,
    help="Upper bound of fixed barycenter support grid",
)
parser.add_argument("--bary_iters", type=int, default=25, help="Sinkhorn barycenter iterations")
parser.add_argument("--outer_epochs", type=int, default=25, help="Maximum outer BCD epochs")
parser.add_argument(
    "--inner_epochs",
    type=int,
    default=3,
    help="Inner epochs for each transport-map block",
)
parser.add_argument(
    "--alpha_epochs",
    type=int,
    default=5,
    help="Alpha-only update epochs after each map-update cycle",
)
parser.add_argument("--hidden_dim", type=int, default=64, help="ICNN hidden width")
parser.add_argument(
    "--nonneg_mode",
    type=str,
    default="relu",
    choices=["relu", "softplus_stable"],
    help="Nonnegative ICNN weight parameterization",
)
parser.add_argument("--lr_map", type=float, default=3e-3, help="Learning rate for ICNN maps")
parser.add_argument(
    "--lr_alpha",
    type=float,
    default=1e-3,
    help="Learning rate for simplex weights",
)
parser.add_argument("--epsilon_bary", type=float, default=0.25, help="Entropic barycenter epsilon")
parser.add_argument(
    "--blur_loss",
    type=float,
    default=0.10,
    help="GeomLoss Sinkhorn blur for prediction objective/evaluation",
)
parser.add_argument("--lambda_penalty", type=float, default=2.0, help="Soft anchor penalty weight")
parser.add_argument(
    "--anchor_mode",
    type=str,
    default="four",
    choices=["two", "four"],
    help="Anchor set: two diagonal anchors or four anchors",
)
parser.add_argument(
    "--response_noise_std",
    type=float,
    default=0.03,
    help="Gaussian perturbation on training response support points",
)
parser.add_argument(
    "--same_cloud",
    action="store_true",
    help="Use the same observed cloud to generate responses; useful for debugging",
)
parser.add_argument(
    "--alpha_true",
    type=str,
    default="0.25,0.35,0.40",
    help="Comma-separated true weights, e.g. 0.25,0.35,0.40",
)
parser.add_argument("--seed", type=int, default=20260603, help="Base random seed")
parser.add_argument("--gpu_id", type=int, default=0, help="GPU id if CUDA is available")
parser.add_argument(
    "--early_stop_tol",
    type=float,
    default=1e-4,
    help="Relative tolerance for validation-loss improvement",
)
parser.add_argument(
    "--early_stop_patience",
    type=int,
    default=4,
    help="Patience for validation-loss early stopping",
)
parser.add_argument(
    "--warmup_epochs",
    type=int,
    default=5,
    help="Minimum epochs before early stopping is allowed",
)
parser.add_argument(
    "--eval_grid_size",
    type=int,
    default=60,
    help="Grid size per side for L2 map-error evaluation",
)
parser.add_argument(
    "--out_csv",
    type=str,
    default="p2_entropic_result.csv",
    help="Raw trial output CSV path",
)
parser.add_argument(
    "--summary_csv",
    type=str,
    default="p2_entropic_summary_all.csv",
    help="Append-only summary output CSV path",
)
args = parser.parse_args()

# Must be set before importing torch.
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from scipy.stats import wishart
from geomloss import SamplesLoss

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def parse_alpha_true(alpha_string: str) -> np.ndarray:
    """Validate and normalize the three population weights (reference first)."""
    alpha = np.asarray([float(x.strip()) for x in alpha_string.split(",")], dtype=np.float32)
    if alpha.shape[0] != 3:
        raise ValueError("For p=2, alpha_true must contain exactly 3 weights: alpha0,alpha1,alpha2.")
    if np.any(alpha < 0):
        raise ValueError("alpha_true must be nonnegative.")
    total = float(alpha.sum())
    if total <= 0:
        raise ValueError("alpha_true must have positive sum.")
    alpha = alpha / total
    return alpha.astype(np.float32)


def set_seed(seed: int) -> None:
    """Set the NumPy and PyTorch random seeds for one replication."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def inv_softplus(y: float) -> float:
    """Convert a positive initialization value to its unconstrained parameter."""
    return math.log(math.expm1(float(y)))


# ------------------------------------------------------------------------------
# Step 2: Learned Convex Potentials and Ground-Truth Transport Maps
# ------------------------------------------------------------------------------

class ICNN(nn.Module):
    """Two-hidden-layer convex potential with nonnegative hidden-path weights."""

    def __init__(
        self,
        dim=2,
        hidden_dim=64,
        nonneg_mode="relu",
        init_positive=0.03,
        identity_strength=0.5,
    ):
        super().__init__()
        self.nonneg_mode = nonneg_mode
        self.identity_strength = identity_strength
        self.W_x0 = nn.Linear(dim, hidden_dim)
        self.W_z1 = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.W_x1 = nn.Linear(dim, hidden_dim)
        self.W_z2 = nn.Linear(hidden_dim, 1, bias=False)
        self.W_x2 = nn.Linear(dim, 1)

        if self.nonneg_mode == "relu":
            nn.init.normal_(self.W_z1.weight, mean=init_positive, std=0.02)
            nn.init.normal_(self.W_z2.weight, mean=init_positive, std=0.02)
        elif self.nonneg_mode == "softplus_stable":
            raw_center = inv_softplus(init_positive)
            nn.init.normal_(self.W_z1.weight, mean=raw_center, std=0.02)
            nn.init.normal_(self.W_z2.weight, mean=raw_center, std=0.02)
        else:
            raise ValueError(f"unknown nonneg_mode={self.nonneg_mode}")

    def positive(self, weight):
        if self.nonneg_mode == "relu":
            return F.relu(weight)
        return F.softplus(weight)

    def forward(self, x):
        z1 = F.softplus(self.W_x0(x))
        z2 = F.softplus(F.linear(z1, self.positive(self.W_z1.weight)) + self.W_x1(x))
        out = F.linear(z2, self.positive(self.W_z2.weight)) + self.W_x2(x)
        return out + self.identity_strength * torch.sum(x ** 2, dim=1, keepdim=True)


def push_forward_icnn(model, x, is_training=True):
    """Evaluate the potential gradient; retain its graph for map updates only."""
    with torch.enable_grad():
        x_in = x.clone().requires_grad_(True)
        phi = model(x_in)
        grad = torch.autograd.grad(
            outputs=phi,
            inputs=x_in,
            grad_outputs=torch.ones_like(phi),
            create_graph=is_training,
            retain_graph=is_training,
        )[0]
    return grad


class GroundTruthMap:
    """Gradient of 0.5*||x||^2 plus a trigonometric perturbation potential.

    The perturbation gradient vanishes at the corners of [-2,2]^2, so the
    ground-truth maps fix the anchor points used below.
    """

    def __init__(self, mode="bubble", strength=0.2):
        self.mode = mode
        self.c = strength

    def push_forward(self, x_tensor):
        with torch.enable_grad():
            x = x_tensor.clone().requires_grad_(True)
            identity_phi = 0.5 * torch.sum(x ** 2, dim=1)
            a = math.pi / 2.0  # matched to [-2, 2]
            if self.mode == "bubble":
                perturb = self.c * torch.cos(a * x[:, 0]) * torch.cos(a * x[:, 1])
            elif self.mode == "wave":
                perturb = self.c * (torch.cos(a * x[:, 0]) + torch.cos(a * x[:, 1]))
            elif self.mode == "cross_wave":
                perturb = self.c * (torch.cos(a * x[:, 0]) - torch.cos(a * x[:, 1]))
            else:
                perturb = 0.0
            phi = identity_phi + perturb
            return torch.autograd.grad(phi.sum(), x)[0]


# ------------------------------------------------------------------------------
# Step 3: Particle Sampling, Numerical Barycenters, and Anchor Points
# ------------------------------------------------------------------------------

def sample_truncated_mvn(mean, cov, num_particles, bounds=(-2.0, 2.0)):
    """Rejection-sample Gaussian particles in a box.

    The original numerical fallback is retained: if no point is accepted
    after 300 batches, uniform particles in the box are used instead.
    """
    pts = []
    lower, upper = bounds
    attempts = 0
    while len(pts) < num_particles:
        attempts += 1
        batch_size = max(num_particles, 256)
        batch = np.random.multivariate_normal(mean, cov, batch_size)
        valid = batch[
            (batch[:, 0] >= lower) & (batch[:, 0] <= upper) &
            (batch[:, 1] >= lower) & (batch[:, 1] <= upper)
        ]
        pts.extend(valid)
        if attempts > 300 and len(pts) == 0:
            pts.extend(np.random.uniform(lower, upper, size=(num_particles, 2)))
    return np.asarray(pts[:num_particles], dtype=np.float32)


def make_grid_support(grid_n, low=-2.4, high=2.4):
    """Construct the fixed numerical response support (grid_n squared points)."""
    grid = torch.linspace(low, high, grid_n, device=device)
    xx, yy = torch.meshgrid(grid, grid, indexing="ij")
    return torch.stack([xx.reshape(-1), yy.reshape(-1)], dim=1)


def entropic_barycenter_weights(support, clouds, weights, epsilon=0.25, n_iters=25):
    """Approximate entropic barycenter masses on a fixed support.

    support: (S, 2)
    clouds: list of K tensors, each (M, 2)
    weights: tensor of length K, nonnegative and summing to one

    Input particles have uniform masses. The fixed number of iterations
    remains differentiable for map/weight updates; it does not compute an
    exact unregularized Wasserstein barycenter.
    """
    tiny = 1e-12
    kernels = []
    masses = []
    scalings_v = []

    for x in clouds:
        cost = torch.cdist(support, x).pow(2)
        kernel = torch.exp(-cost / epsilon).clamp_min(tiny)
        kernels.append(kernel)
        masses.append(torch.full((x.shape[0],), 1.0 / x.shape[0], device=x.device, dtype=x.dtype))
        scalings_v.append(torch.ones(x.shape[0], device=x.device, dtype=x.dtype))

    q = torch.full(
        (support.shape[0],),
        1.0 / support.shape[0],
        device=support.device,
        dtype=support.dtype,
    )
    for _ in range(n_iters):
        logs = []
        new_v = []
        for kernel, mass, v in zip(kernels, masses, scalings_v):
            u = q / (kernel @ v).clamp_min(tiny)
            v = mass / (kernel.T @ u).clamp_min(tiny)
            new_v.append(v)
            logs.append(torch.log((kernel @ v).clamp_min(tiny)))
        scalings_v = new_v
        log_q = sum(weights[j] * logs[j] for j in range(len(clouds)))
        q = torch.softmax(log_q, dim=0)

    return q


def make_gt_maps():
    """Return T0*, T1*, T2* with their simulation perturbation strengths."""
    return [
        GroundTruthMap("bubble", 0.20),
        GroundTruthMap("wave", -0.16),
        GroundTruthMap("cross_wave", 0.14),
    ]


def anchor_targets(gt_maps, anchor_mode="four"):
    """Construct four corner anchors or two opposite-corner anchors and targets."""
    if anchor_mode == "four":
        points = [[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]]
    else:
        points = [[-2.0, -2.0], [2.0, 2.0]]
    anchors_src = torch.tensor(points, dtype=torch.float32, device=device)
    with torch.no_grad():
        anchors_tgt = [gt.push_forward(anchors_src).detach() for gt in gt_maps]
    return anchors_src, anchors_tgt


# ------------------------------------------------------------------------------
# Step 4: Generate Training, Validation, and Test Distributions
# ------------------------------------------------------------------------------

def generate_p2_dataset(
    num_subjects,
    m_particles,
    gt_maps,
    alpha_true,
    support,
    noise_std,
    is_train=True,
    same_cloud=False,
):
    """Generate p=2 data: one reference component plus two predictor components.

    Each subject has independently generated predictor parameters. Within a
    subject, the two predictors also have independent parameter draws.
    Observed inputs and hidden clouds share those parameters, but use separate
    particle draws unless same_cloud=True (a debugging option).

    Clean responses are finite-iteration numerical barycenters. Training
    perturbs their support points while keeping their masses fixed; validation
    and test data use the unperturbed support (is_train=False).
    """
    data = []
    bounds = (-2.0, 2.0)
    alpha = torch.tensor(alpha_true, dtype=torch.float32, device=device)

    for _ in range(num_subjects):
        x_obs = []
        x_hidden = []

        for j in range(3):
            if j == 0:
                mean = np.zeros(2)
                cov = np.eye(2)
            else:
                mean = np.random.uniform(-1.0, 1.0, size=2)
                # SciPy uses a scale matrix: E[cov] = df * scale = 1.2 * I2.
                cov = wishart.rvs(df=4, scale=np.eye(2) * 0.3)

            obs = torch.tensor(sample_truncated_mvn(mean, cov, m_particles, bounds), device=device)
            x_obs.append(obs)
            if same_cloud:
                x_hidden.append(obs.clone())
            else:
                x_hidden.append(torch.tensor(sample_truncated_mvn(mean, cov, m_particles, bounds), device=device))

        with torch.no_grad():
            pushed = [gt_maps[j].push_forward(x_hidden[j]) for j in range(3)]
            q_clean = entropic_barycenter_weights(
                support,
                pushed,
                alpha,
                epsilon=args.epsilon_bary,
                n_iters=args.bary_iters,
            ).detach()

            actual_noise = noise_std if is_train else 0.0
            y_support = support + actual_noise * torch.randn_like(support)

        data.append({
            "X": x_obs,
            "q_obs": q_clean.detach(),
            "Y_support": y_support.detach(),
            "is_train": is_train,
        })

    return data


# ------------------------------------------------------------------------------
# Step 5: Prediction and Parameter Diagnostics
# ------------------------------------------------------------------------------

def evaluate_l2_error_map(model, gt_map, bounds=(-2.0, 2.0), grid_size=60):
    """Return the root mean squared Euclidean map error on a regular grid."""
    model.eval()
    lo, hi = bounds
    grid_val = np.linspace(lo, hi, grid_size)
    xx, yy = np.meshgrid(grid_val, grid_val)
    x_grid = torch.tensor(np.c_[xx.ravel(), yy.ravel()], dtype=torch.float32, device=device)

    y_gt = gt_map.push_forward(x_grid).detach().cpu().numpy()
    y_pred = push_forward_icnn(model, x_grid, is_training=False).detach().cpu().numpy()

    return float(np.sqrt(np.mean(np.linalg.norm(y_pred - y_gt, axis=1) ** 2)))


def evaluate_prediction_loss(
    models, data, alpha_logits, support, loss_fn, epsilon_bary, bary_iters,
):
    """Average Sinkhorn prediction divergence across subjects, without a root."""
    for m in models:
        m.eval()
    alpha = torch.softmax(alpha_logits, dim=0)
    total = 0.0

    for batch in data:
        clouds = []
        for j in range(3):
            clouds.append(push_forward_icnn(models[j], batch["X"][j], is_training=False).detach())
        q_pred = entropic_barycenter_weights(
            support,
            clouds,
            alpha,
            epsilon=epsilon_bary,
            n_iters=bary_iters,
        )
        total += float(loss_fn(q_pred, support, batch["q_obs"], batch["Y_support"]).detach().cpu())

    return total / max(len(data), 1)


def evaluate_global_anchor_penalty(models, anchors_src, anchors_tgt):
    """Average squared anchor discrepancy across points and the three maps."""
    vals = []
    for j, model in enumerate(models):
        model.eval()
        pred = push_forward_icnn(model, anchors_src, is_training=False).detach()
        vals.append(torch.mean(torch.sum((pred - anchors_tgt[j]) ** 2, dim=1)))
    return float(torch.stack(vals).mean().detach().cpu())


# ------------------------------------------------------------------------------
# Step 6: Validation-Based Stopping and Model Checkpoints
# ------------------------------------------------------------------------------

class EarlyStopper:
    """Track validation improvements scaled by max(1, abs(best)) after warm-up."""

    def __init__(self, patience=4, min_rel_improve=1e-4, warmup=5):
        self.patience = patience
        self.min_rel_improve = min_rel_improve
        self.warmup = warmup
        self.best = float("inf")
        self.best_epoch = -1
        self.counter = 0

    def step(self, value, epoch):
        if value < self.best:
            rel_improve = (self.best - value) / max(1.0, abs(self.best))
        else:
            rel_improve = -float("inf")

        improved = (epoch < self.warmup) or (rel_improve > self.min_rel_improve)
        if improved:
            if value < self.best:
                self.best = value
                self.best_epoch = epoch
            self.counter = 0
            return False

        self.counter += 1
        return self.counter >= self.patience


def clone_state(models, alpha_logits):
    """Copy all map parameters and simplex logits for later evaluation."""
    return {
        "models": [copy.deepcopy(m.state_dict()) for m in models],
        "alpha_logits": alpha_logits.detach().clone(),
    }


def restore_state(models, alpha_logits, state):
    """Restore a saved validation-selected state without changing optimizers."""
    for m, sd in zip(models, state["models"]):
        m.load_state_dict(sd)
    with torch.no_grad():
        alpha_logits.copy_(state["alpha_logits"])


# ------------------------------------------------------------------------------
# Step 7: Fit MTDR and Evaluate One Monte Carlo Replication
# ------------------------------------------------------------------------------

def train_one_trial(trial_id: int, seed: int):
    """Generate data, alternate map/weight blocks, and evaluate a saved fit."""
    set_seed(seed)
    alpha_true = parse_alpha_true(args.alpha_true)
    support = make_grid_support(args.grid_n, args.support_low, args.support_high)
    gt_maps = make_gt_maps()
    anchors_src, anchors_tgt = anchor_targets(gt_maps, args.anchor_mode)

    print(
        f"=== p=2 entropic MTDR | trial={trial_id} | seed={seed} | device={device} | "
        f"N={args.N}, N_val={args.N_val}, N_test={args.N_test}, M={args.M}, grid_n={args.grid_n}, "
        f"outer={args.outer_epochs}, inner={args.inner_epochs}, alpha_epochs={args.alpha_epochs}, "
        f"anchor_mode={args.anchor_mode}, lambda={args.lambda_penalty} ==="
    )

    data_start = time.time()
    train_data = generate_p2_dataset(
        args.N,
        args.M,
        gt_maps,
        alpha_true,
        support,
        args.response_noise_std,
        is_train=True,
        same_cloud=args.same_cloud,
    )
    val_data = generate_p2_dataset(
        args.N_val,
        args.M,
        gt_maps,
        alpha_true,
        support,
        args.response_noise_std,
        is_train=False,
        same_cloud=args.same_cloud,
    )
    test_data = generate_p2_dataset(
        args.N_test,
        args.M,
        gt_maps,
        alpha_true,
        support,
        args.response_noise_std,
        is_train=False,
        same_cloud=args.same_cloud,
    )
    data_seconds = time.time() - data_start
    print(f">>> Data generated in {data_seconds:.2f}s")

    models = [ICNN(hidden_dim=args.hidden_dim, nonneg_mode=args.nonneg_mode).to(device) for _ in range(3)]
    alpha_logits = nn.Parameter(torch.zeros(3, device=device))
    # Zero logits initialize all three estimated weights to 1/3.

    opts = [optim.Adam(models[j].parameters(), lr=args.lr_map) for j in range(3)]
    opt_alpha = optim.Adam([alpha_logits], lr=args.lr_alpha)
    loss_fn = SamplesLoss(loss="sinkhorn", p=2, blur=args.blur_loss)

    early = EarlyStopper(
        patience=args.early_stop_patience,
        min_rel_improve=args.early_stop_tol,
        warmup=args.warmup_epochs,
    )
    best_state = clone_state(models, alpha_logits)

    last_train_total = float("nan")
    last_train_pred = float("nan")
    last_train_anchor = float("nan")
    last_alpha_phase_loss = float("nan")
    last_val_loss = float("nan")
    last_global_anchor = float("nan")
    epochs_run = 0
    train_start = time.time()

    for epoch in range(args.outer_epochs):
        epochs_run = epoch + 1
        epoch_total_sum = 0.0
        epoch_pred_sum = 0.0
        epoch_anchor_sum = 0.0
        epoch_count = 0

        # ---------------------------------------------------------
        # Map blocks: update one ICNN at a time, alpha fixed.
        # ---------------------------------------------------------
        for block in range(3):
            for j in range(3):
                models[j].train(j == block)

            for _ in range(args.inner_epochs):
                np.random.shuffle(train_data)
                for batch in train_data:
                    opts[block].zero_grad(set_to_none=True)

                    clouds = []
                    for j in range(3):
                        if j == block:
                            clouds.append(push_forward_icnn(models[j], batch["X"][j], is_training=True))
                        else:
                            clouds.append(push_forward_icnn(models[j], batch["X"][j], is_training=False).detach())

                    alpha_fixed = torch.softmax(alpha_logits, dim=0).detach()
                    q_pred = entropic_barycenter_weights(
                        support,
                        clouds,
                        alpha_fixed,
                        epsilon=args.epsilon_bary,
                        n_iters=args.bary_iters,
                    )

                    loss_pred = loss_fn(q_pred, support, batch["q_obs"], batch["Y_support"])

                    if args.lambda_penalty > 0:
                        t_anchor = push_forward_icnn(models[block], anchors_src, is_training=True)
                        loss_anchor = torch.mean(torch.sum((t_anchor - anchors_tgt[block]) ** 2, dim=1))
                    else:
                        loss_anchor = torch.tensor(0.0, device=device)

                    loss = loss_pred + args.lambda_penalty * loss_anchor
                    loss.backward()
                    opts[block].step()

                    epoch_total_sum += float(loss.detach().cpu())
                    epoch_pred_sum += float(loss_pred.detach().cpu())
                    epoch_anchor_sum += float(loss_anchor.detach().cpu())
                    epoch_count += 1

        # ---------------------------------------------------------
        # Alpha block: update simplex weights only, maps fixed.
        # ---------------------------------------------------------
        alpha_loss_sum = 0.0
        alpha_loss_count = 0
        for _ in range(args.alpha_epochs):
            np.random.shuffle(train_data)
            for batch in train_data:
                opt_alpha.zero_grad(set_to_none=True)
                clouds = []
                for j in range(3):
                    clouds.append(push_forward_icnn(models[j], batch["X"][j], is_training=False).detach())

                alpha = torch.softmax(alpha_logits, dim=0)
                q_pred = entropic_barycenter_weights(
                    support,
                    clouds,
                    alpha,
                    epsilon=args.epsilon_bary,
                    n_iters=args.bary_iters,
                )

                loss_alpha = loss_fn(q_pred, support, batch["q_obs"], batch["Y_support"])
                loss_alpha.backward()
                opt_alpha.step()

                alpha_loss_sum += float(loss_alpha.detach().cpu())
                alpha_loss_count += 1

        last_train_total = epoch_total_sum / max(epoch_count, 1)
        last_train_pred = epoch_pred_sum / max(epoch_count, 1)
        last_train_anchor = epoch_anchor_sum / max(epoch_count, 1)
        last_alpha_phase_loss = alpha_loss_sum / max(alpha_loss_count, 1)
        last_val_loss = evaluate_prediction_loss(
            models,
            val_data,
            alpha_logits,
            support,
            loss_fn,
            args.epsilon_bary,
            args.bary_iters,
        )
        last_global_anchor = evaluate_global_anchor_penalty(models, anchors_src, anchors_tgt)
        alpha_now = torch.softmax(alpha_logits, dim=0).detach().cpu().numpy()

        print(
            f"Epoch {epochs_run:02d}/{args.outer_epochs} | "
            f"train_total={last_train_total:.5f} | train_pred={last_train_pred:.5f} | "
            f"train_anchor={last_train_anchor:.5f} | alpha_phase={last_alpha_phase_loss:.5f} | "
            f"val_pred={last_val_loss:.5f} | global_anchor={last_global_anchor:.5f} | "
            f"alpha={np.round(alpha_now, 4)}"
        )

        if last_val_loss < early.best:
            best_state = clone_state(models, alpha_logits)

        if early.step(last_val_loss, epoch):
            print(f">>> Early stopping triggered at epoch {epochs_run}; best epoch={early.best_epoch + 1}")
            break

    restore_state(models, alpha_logits, best_state)
    # Final test evaluation uses the saved state, not the last optimization step.
    train_seconds = time.time() - train_start

    val_loss = evaluate_prediction_loss(
        models,
        val_data,
        alpha_logits,
        support,
        loss_fn,
        args.epsilon_bary,
        args.bary_iters,
    )
    test_loss = evaluate_prediction_loss(
        models,
        test_data,
        alpha_logits,
        support,
        loss_fn,
        args.epsilon_bary,
        args.bary_iters,
    )
    global_anchor = evaluate_global_anchor_penalty(models, anchors_src, anchors_tgt)

    alpha_hat = torch.softmax(alpha_logits, dim=0).detach().cpu().numpy()
    alpha_l1_err = float(np.sum(np.abs(alpha_hat - alpha_true)))
    alpha_l2_err = float(np.sqrt(np.sum((alpha_hat - alpha_true) ** 2)))

    map_errs = [
        evaluate_l2_error_map(models[j], gt_maps[j], grid_size=args.eval_grid_size)
        for j in range(3)
    ]

    result = {
        "trial_id": trial_id,
        "seed": seed,
        "N": args.N,
        "N_val": args.N_val,
        "N_test": args.N_test,
        "M": args.M,
        "grid_n": args.grid_n,
        "support_size": args.grid_n ** 2,
        "bary_iters": args.bary_iters,
        "outer_epochs": args.outer_epochs,
        "inner_epochs": args.inner_epochs,
        "alpha_epochs": args.alpha_epochs,
        "epochs_run": epochs_run,
        "best_epoch": early.best_epoch + 1,
        "hidden_dim": args.hidden_dim,
        "nonneg_mode": args.nonneg_mode,
        "epsilon_bary": args.epsilon_bary,
        "blur_loss": args.blur_loss,
        "lambda_penalty": args.lambda_penalty,
        "anchor_mode": args.anchor_mode,
        "response_noise_std_train": args.response_noise_std,
        "test_response_clean": True,
        "same_cloud": args.same_cloud,
        "alpha0_true": alpha_true[0],
        "alpha1_true": alpha_true[1],
        "alpha2_true": alpha_true[2],
        "alpha0_hat": alpha_hat[0],
        "alpha1_hat": alpha_hat[1],
        "alpha2_hat": alpha_hat[2],
        "alpha_l1_err": alpha_l1_err,
        "alpha_l2_err": alpha_l2_err,
        "T0_L2": map_errs[0],
        "T1_L2": map_errs[1],
        "T2_L2": map_errs[2],
        "last_train_total_loss": last_train_total,
        "last_train_pred_loss": last_train_pred,
        "last_train_anchor_loss": last_train_anchor,
        "last_alpha_phase_loss": last_alpha_phase_loss,
        "val_sinkhorn_loss": val_loss,
        "test_sinkhorn_loss": test_loss,
        "global_anchor_penalty": global_anchor,
        "data_seconds": data_seconds,
        "train_seconds": train_seconds,
        "total_seconds": data_seconds + train_seconds,
    }

    print("\n>>> Trial result")
    print(pd.DataFrame([result]).to_string(index=False))
    return result


# ------------------------------------------------------------------------------
# Step 8: Monte Carlo Replications and CSV Summaries
# ------------------------------------------------------------------------------

def run_trials():
    """Run replications sequentially and summarize successful trials only."""
    results = []
    start = time.time()

    for t in range(args.trials):
        trial_seed = args.seed + 1000 * t
        try:
            res = train_one_trial(t, trial_seed)
            results.append(res)
            pd.DataFrame(results).to_csv(args.out_csv, index=False)
            print(f">>> Saved running results to {args.out_csv}")
        except Exception as exc:
            print(f"❌ Trial {t} failed: {repr(exc)}")
            if args.trials == 1:
                raise
        finally:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    if not results:
        raise RuntimeError("No successful trials. Please check errors above.")

    df = pd.DataFrame(results)
    summary = {
        "trials_successful": len(df),
        "N": args.N,
        "N_val": args.N_val,
        "N_test": args.N_test,
        "M": args.M,
        "alpha_true": args.alpha_true,
        "alpha0_true": df["alpha0_true"].iloc[0],
        "alpha1_true": df["alpha1_true"].iloc[0],
        "alpha2_true": df["alpha2_true"].iloc[0],
        "grid_n": args.grid_n,
        "support_size": args.grid_n ** 2,
        "support_low": args.support_low,
        "support_high": args.support_high,
        "bary_iters": args.bary_iters,
        "epsilon_bary": args.epsilon_bary,
        "blur_loss": args.blur_loss,
        "outer_epochs": args.outer_epochs,
        "inner_epochs": args.inner_epochs,
        "alpha_epochs": args.alpha_epochs,
        "hidden_dim": args.hidden_dim,
        "nonneg_mode": args.nonneg_mode,
        "lr_map": args.lr_map,
        "lr_alpha": args.lr_alpha,
        "lambda_penalty": args.lambda_penalty,
        "anchor_mode": args.anchor_mode,
        "response_noise_std": args.response_noise_std,
        "same_cloud": args.same_cloud,
        "alpha_l1_err_mean": df["alpha_l1_err"].mean(),
        "alpha_l1_err_std": df["alpha_l1_err"].std(),
        "alpha_l2_err_mean": df["alpha_l2_err"].mean(),
        "alpha_l2_err_std": df["alpha_l2_err"].std(),
        "T0_L2_mean": df["T0_L2"].mean(),
        "T0_L2_std": df["T0_L2"].std(),
        "T1_L2_mean": df["T1_L2"].mean(),
        "T1_L2_std": df["T1_L2"].std(),
        "T2_L2_mean": df["T2_L2"].mean(),
        "T2_L2_std": df["T2_L2"].std(),
        "val_sinkhorn_loss_mean": df["val_sinkhorn_loss"].mean(),
        "val_sinkhorn_loss_std": df["val_sinkhorn_loss"].std(),
        "test_sinkhorn_loss_mean": df["test_sinkhorn_loss"].mean(),
        "test_sinkhorn_loss_std": df["test_sinkhorn_loss"].std(),
        "global_anchor_penalty_mean": df["global_anchor_penalty"].mean(),
        "epochs_run_mean": df["epochs_run"].mean(),
        "total_minutes": (time.time() - start) / 60.0,
    }
    summary_df = pd.DataFrame([summary])
    # Preserve existing summary columns when appending another configuration.
    if os.path.exists(args.summary_csv):
        existing = pd.read_csv(args.summary_csv)
        for col in summary_df.columns:
            if col not in existing.columns:
                existing[col] = np.nan
        for col in existing.columns:
            if col not in summary_df.columns:
                summary_df[col] = np.nan
        summary_df = summary_df[existing.columns]
        pd.concat([existing, summary_df], ignore_index=True).to_csv(args.summary_csv, index=False)
    else:
        summary_df.to_csv(args.summary_csv, index=False)

    print("\n=== All trials completed ===")
    print(summary_df.to_string(index=False))
    print(f"\n>>> Raw results saved to {args.out_csv}")
    print(f">>> Summary appended to {args.summary_csv}")


if __name__ == "__main__":
    run_trials()
