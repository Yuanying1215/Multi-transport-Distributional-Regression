# ==============================================================================
# Script: Table1.py
# Description: Simulation study for the Single Predictor setting (Section 6.1).
#              Fits MTDR to two-dimensional distributions using input-convex
#              neural networks (ICNNs). One transported reference distribution
#              and one transported random predictor are aggregated through
#              optimal matching and displacement interpolation.
#
# Method: Multi-transport Distributional Regression (MTDR)
# Key Algorithm: Block-coordinate optimization of T0, T1, and the reference weight.
#
# Output:
#   1. mtdr_icnn_alpha{alpha}_N{N}_M{M}.csv: Trial-level results for one setting.
#   2. summary_all.csv: Monte Carlo means and sample SDs, appended across runs.
#      The Test_W2sq and Test_RMSE summaries support the prediction comparison
#      in Table 1. The script does not assemble or round the manuscript table.
#
# Data:
#   Synthetic data are generated within each replication; no data files needed.
#   The reference is a truncated N(0, I2) law. Each predictor has a random mean
#   and Wishart covariance, with observed and hidden particle clouds sampled
#   independently by default. Training responses have Gaussian particle noise;
#   validation and test responses are noise-free numerical targets.
#
# Dependencies:
#   Python 3, numpy, pandas, scipy, torch (PyTorch), geomloss.
#   CUDA is used when available; otherwise the script uses the CPU.
#
# Usage (run from this folder; use --help to list all options):
#   Quick execution check, NOT a full Monte Carlo experiment:
#     python Table1.py --alpha 0.5 --N 20 --M 100 --trials 2 \
#         --max_workers 1 --outer_epochs 10 --inner_epochs 5 --alpha_epochs 2
#
#   Example configuration on one GPU:
#     python Table1.py --alpha 0.5 --N 50 --M 400 --trials 50 \
#         --gpu_id 0 --max_workers 1 --outdir results_alpha05_N50_M400
#
# Parameter Conventions:
#   --alpha is alpha0*, the REFERENCE weight; alpha1* = 1 - alpha0*.
#   --N is the training sample size n; --M is the particles per distribution m.
#   Validation and test sizes are max(2, round(N * ratio)), with default ratios
#   0.2 and 0.3. Replication t uses seed t * seed_stride (default stride: 1000).
#
# Evaluation and Numerical Notes:
#   - Training uses GeomLoss Sinkhorn divergence (p=2, default blur=0.05).
#   - Validation and testing use exact balanced empirical W2^2, computed by
#     Hungarian matching between equal-size, equal-weight particle clouds.
#   - Test_W2sq is the test mean squared distance; Test_RMSE is its square root.
#     Test_RMSE_Mean averages replication-level RMSEs, not squared losses.
#   - Err_T0/Err_T1 are grid-based root-mean-squared map discrepancies, not
#     unnormalized Lebesgue L2 integrals. Weight/map errors are diagnostics;
#     an inactive map is not identifiable when its true weight is zero.
#   - [-2,2]^2 bounds the sampled inputs. Gaussian response noise and learned
#     ICNN maps are not clipped; soft corner penalties do not enforce bounds.
#   - Early stopping uses validation data; evaluation restores the saved best
#     validation state. The test data are not used for model selection.
#   - Successful trials are summarized; inspect Trials_Success and error logs.
#     Repeating a setting overwrites its raw CSV and appends another summary.
#     Use a separate --outdir for each intended experiment/run.
#   - Seeds are fixed per trial, but bitwise reproducibility across devices and
#     software versions is not guaranteed. This script is intended for CLI use:
#     argument parsing and device setup occur at module level.
# ==============================================================================

# ------------------------------------------------------------------------------
# Step 1: Runtime Configuration and Command-Line Options
# ------------------------------------------------------------------------------
# Configure numerical-library threads and select the GPU before importing torch.

import os
import sys

# Limit BLAS/OpenMP threads before importing numpy/scipy/torch.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Avoid a script named code.py shadowing Python's stdlib code module.
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path = [p for p in sys.path if os.path.abspath(p or os.getcwd()) != _script_dir]

import argparse

parser = argparse.ArgumentParser(description="MTDR p=1 2D ICNN Monte Carlo Simulation")
parser.add_argument(
    "--alpha",
    type=float,
    required=True,
    help="True alpha*, e.g. 0, 0.25, 0.5, 0.75, 1",
)
parser.add_argument("--N", type=int, required=True, help="Number of training distributions")
parser.add_argument("--trials", type=int, default=50, help="Number of Monte Carlo trials")
parser.add_argument("--M", type=int, default=400, help="Particles per distribution")
parser.add_argument("--gpu_id", type=int, default=0, help="GPU id to use")
parser.add_argument(
    "--max_workers",
    type=int,
    default=None,
    help="Parallel MC workers. Default: 1 on CUDA, up to 8 on CPU",
)

# Optimization
parser.add_argument("--outer_epochs", type=int, default=40, help="Maximum BCD outer epochs")
parser.add_argument("--inner_epochs", type=int, default=20, help="Inner epochs for each map block")
parser.add_argument("--alpha_epochs", type=int, default=5, help="Epochs for alpha-only block")
parser.add_argument("--lr_map", type=float, default=3e-3, help="Learning rate for ICNN maps")
parser.add_argument("--lr_alpha", type=float, default=1e-3, help="Learning rate for alpha logit")
parser.add_argument("--blur_train", type=float, default=0.05, help="Sinkhorn blur for training")
parser.add_argument("--lambda_penalty", type=float, default=2.0, help="Anchor soft-penalty weight")
parser.add_argument(
    "--anchor_mode",
    choices=["four", "three"],
    default="four",
    help="four: boundry anchors; three: simplex anchors",
)

# Early stopping
parser.add_argument(
    "--val_ratio",
    type=float,
    default=0.2,
    help="Validation set size as a fraction of N",
)
parser.add_argument(
    "--test_ratio",
    type=float,
    default=0.3,
    help="Test set size as a fraction of N",
)
parser.add_argument("--early_stop_patience", type=int, default=5, help="Validation patience")
parser.add_argument(
    "--early_stop_tol",
    type=float,
    default=1e-3,
    help="Minimum relative validation improvement",
)
parser.add_argument(
    "--warmup_epochs",
    type=int,
    default=5,
    help="No early stop before this many outer epochs",
)
parser.add_argument(
    "--val_every",
    type=int,
    default=1,
    help="Evaluate validation every k outer epochs",
)

# Data
parser.add_argument("--noise_std", type=float, default=0.05, help="Training response noise std")
parser.add_argument("--bounds_low", type=float, default=-2.0, help="Lower bound of compact domain")
parser.add_argument("--bounds_high", type=float, default=2.0, help="Upper bound of compact domain")
parser.add_argument(
    "--hidden_independent",
    action="store_true",
    default=True,
    help="Use independent hidden clouds for response generation",
)
parser.add_argument(
    "--same_cloud_debug",
    action="store_true",
    help="Use observed clouds as hidden clouds; useful for debugging only",
)
parser.add_argument(
    "--eval_grid",
    type=int,
    default=60,
    help="Grid size per dimension for L2 map error",
)
parser.add_argument("--hidden_dim", type=int, default=128, help="ICNN hidden dimension")
parser.add_argument("--outdir", type=str, default=".", help="Directory to write CSV files")
parser.add_argument("--seed_stride", type=int, default=1000, help="Seed offset between trials")
args = parser.parse_args()

# Must be set before importing torch.
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

import math
import time
import copy
import multiprocessing as mp
import concurrent.futures

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from scipy.stats import wishart
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from geomloss import SamplesLoss


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BOUNDS = (float(args.bounds_low), float(args.bounds_high))
if not (abs(BOUNDS[0] + 2.0) < 1e-12 and abs(BOUNDS[1] - 2.0) < 1e-12):
    print(f"[Warning] This script is tuned for [-2, 2]. Current bounds={BOUNDS}.")

print(
    f"=== MTDR ICNN MC | alpha={args.alpha} | N={args.N} | M={args.M} | "
    f"trials={args.trials} | device={device} | bounds={BOUNDS} ===",
    flush=True,
)


# ------------------------------------------------------------------------------
# Step 2: ICNN Potentials and Transport Maps
# ------------------------------------------------------------------------------
# Each learned map is the gradient of a convex potential. The hidden-to-hidden
# and hidden-to-output weights are nonnegative; affine input skip terms are free.
class ICNN(nn.Module):
    """Two hidden layers with Softplus activations and a quadratic potential term."""

    def __init__(self, dim=2, hidden_dim=128):
        super().__init__()
        self.W_x0 = nn.Linear(dim, hidden_dim)
        self.W_z1 = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.W_x1 = nn.Linear(dim, hidden_dim)
        self.W_z2 = nn.Linear(hidden_dim, 1, bias=False)
        self.W_x2 = nn.Linear(dim, 1)

    def forward(self, x):
        z1 = F.softplus(self.W_x0(x))
        # Enforce nonnegative hidden weights using the original ReLU parameterization.
        W_z1_pos = F.relu(self.W_z1.weight)
        z2 = F.softplus(F.linear(z1, W_z1_pos) + self.W_x1(x))
        W_z2_pos = F.relu(self.W_z2.weight)
        out = F.linear(z2, W_z2_pos) + self.W_x2(x)
        strong_convex_term = 0.01 * torch.sum(x ** 2, dim=1, keepdim=True)
        return out + strong_convex_term


def push_forward_icnn_constrained(icnn_model, x, is_training=True):
    """Evaluate grad(phi); retain its derivative graph only for map optimization.

    Despite the legacy function name, there is no hard support projection.
    The landmark penalty is added separately in the training loop.
    """
    with torch.enable_grad():
        x_in = x.clone().requires_grad_(True)
        phi_x = icnn_model(x_in)
        grad_x = torch.autograd.grad(
            outputs=phi_x,
            inputs=x_in,
            grad_outputs=torch.ones_like(phi_x),
            create_graph=is_training,
            retain_graph=is_training,
        )[0]
    return grad_x


class GroundTruthMap:
    """
    Ground-truth optimal transport map given by a convex potential gradient.

    At the default L=2, kappa=pi/2 makes the perturbation gradient vanish at
    the four corner landmarks. This is not a claim that every boundary point
    is fixed. The bubble/wave coefficients are specified in train_one_trial().
    """
    def __init__(self, mode="bubble", strength=0.2, L=2.0):
        self.mode = mode
        self.c = float(strength)
        self.L = float(L)
        self.kappa = math.pi / self.L

    def push_forward_constrained(self, X_tensor):
        with torch.enable_grad():
            X = X_tensor.clone().requires_grad_(True)
            identity_phi = 0.5 * torch.sum(X ** 2, dim=1)
            k = self.kappa

            if self.mode == "bubble":
                perturbation = self.c * torch.cos(k * X[:, 0]) * torch.cos(k * X[:, 1])
            elif self.mode == "wave":
                perturbation = self.c * torch.cos(k * X[:, 0]) + self.c * torch.cos(k * X[:, 1])
            else:
                perturbation = torch.zeros_like(identity_phi)

            phi = identity_phi + perturbation
            T_x = torch.autograd.grad(phi.sum(), X)[0]
        return T_x


# ------------------------------------------------------------------------------
# Step 3: Predictor, Reference, and Response Generation
# ------------------------------------------------------------------------------
def sample_truncated_mvn(mean, cov, num_particles, bounds=BOUNDS):
    """Sample a bivariate normal law truncated to the square by rejection."""
    pts = []
    lower, upper = bounds
    # Rejection sampling. Batch size adapts to avoid many tiny loops.
    batch_size = max(num_particles, 200)
    while len(pts) < num_particles:
        batch = np.random.multivariate_normal(mean, cov, batch_size)
        valid = batch[
            (batch[:, 0] >= lower) & (batch[:, 0] <= upper) &
            (batch[:, 1] >= lower) & (batch[:, 1] <= upper)
        ]
        if valid.shape[0] > 0:
            pts.extend(valid.tolist())
    return np.asarray(pts[:num_particles], dtype=np.float32)


def make_single_sample(
    M_particles, gt_map0, gt_map1, alpha_true, noisy_response, noise_std,
    bounds=BOUNDS, same_cloud_debug=False,
):
    """Generate one pair of input clouds and its clean/noisy response cloud."""
    # A common reference law is represented by independently sampled clouds
    # for each subject; it is not a subject-specific random reference law.
    mean_0 = np.zeros(2)
    cov_0 = np.eye(2)

    X0_obs_np = sample_truncated_mvn(mean_0, cov_0, M_particles, bounds)
    X0_obs = torch.tensor(X0_obs_np, dtype=torch.float32, device=device)
    if same_cloud_debug:
        X0_hidden = X0_obs.clone()
    else:
        X0_hidden_np = sample_truncated_mvn(mean_0, cov_0, M_particles, bounds)
        X0_hidden = torch.tensor(X0_hidden_np, dtype=torch.float32, device=device)

    mean_1 = np.random.uniform(-1.0, 1.0, size=2)
    # SciPy uses the scale-matrix convention: E[cov_1] = 4 * 0.3 * I2.
    cov_1 = wishart.rvs(df=4, scale=np.eye(2) * 0.3)

    X1_obs_np = sample_truncated_mvn(mean_1, cov_1, M_particles, bounds)
    X1_obs = torch.tensor(X1_obs_np, dtype=torch.float32, device=device)
    if same_cloud_debug:
        X1_hidden = X1_obs.clone()
    else:
        X1_hidden_np = sample_truncated_mvn(mean_1, cov_1, M_particles, bounds)
        X1_hidden = torch.tensor(X1_hidden_np, dtype=torch.float32, device=device)

    Y0_hat = gt_map0.push_forward_constrained(X0_hidden).detach()
    Y1_hat = gt_map1.push_forward_constrained(X1_hidden).detach()

    C = cdist(Y0_hat.cpu().numpy(), Y1_hat.cpu().numpy(), metric="sqeuclidean")
    _, col_ind = linear_sum_assignment(C)
    Y1_matched = Y1_hat[col_ind]

    # Optimal pairing of the two empirical measures gives their weighted
    # Wasserstein Frechet mean by displacement interpolation.
    Y_clean = alpha_true * Y0_hat + (1.0 - alpha_true) * Y1_matched
    if noisy_response and noise_std > 0:
        Y_obs = Y_clean + torch.randn_like(Y_clean) * noise_std
    else:
        Y_obs = Y_clean.clone()

    return {"X0": X0_obs, "X1": X1_obs, "Y_obs": Y_obs.detach(), "Y_clean": Y_clean.detach()}


def generate_dataset_splits(
    N_train, N_val, N_test, M_particles, gt_map0, gt_map1, alpha_true,
    noise_std=0.05, bounds=BOUNDS, same_cloud_debug=False,
):
    """Generate independent training, validation, and test subjects."""
    train_data = [
        make_single_sample(
            M_particles,
            gt_map0,
            gt_map1,
            alpha_true,
            True,
            noise_std,
            bounds,
            same_cloud_debug,
        )
        for _ in range(N_train)
    ]
    # Validation/test are clean in this simulation so early stopping monitors the target prediction error.
    val_data = [
        make_single_sample(
            M_particles,
            gt_map0,
            gt_map1,
            alpha_true,
            False,
            0.0,
            bounds,
            same_cloud_debug,
        )
        for _ in range(N_val)
    ]
    test_data = [
        make_single_sample(
            M_particles,
            gt_map0,
            gt_map1,
            alpha_true,
            False,
            0.0,
            bounds,
            same_cloud_debug,
        )
        for _ in range(N_test)
    ]
    return train_data, val_data, test_data


# ------------------------------------------------------------------------------
# Step 4: Prediction Metrics and Map-Error Diagnostics
# ------------------------------------------------------------------------------
def exact_w2_sq_equal_weights(X, Y):
    """Exact balanced W2^2 for two equal-size equal-weight point clouds."""
    X_np = X.detach().cpu().numpy()
    Y_np = Y.detach().cpu().numpy()
    C = cdist(X_np, Y_np, metric="sqeuclidean")
    row_ind, col_ind = linear_sum_assignment(C)
    return float(C[row_ind, col_ind].mean())


def make_predicted_barycenter(model0, model1, X0, X1, alpha_value):
    """Transport both observed clouds, optimally pair them, and interpolate."""
    X0_hat = push_forward_icnn_constrained(model0, X0, is_training=False).detach()
    X1_hat = push_forward_icnn_constrained(model1, X1, is_training=False).detach()
    C = cdist(X0_hat.cpu().numpy(), X1_hat.cpu().numpy(), metric="sqeuclidean")
    _, col_ind = linear_sum_assignment(C)
    Y_pred = alpha_value * X0_hat + (1.0 - alpha_value) * X1_hat[col_ind]
    return Y_pred.detach()


def evaluate_prediction_exact(model0, model1, data, alpha_value):
    """Return test/validation mean empirical W2^2 and its square root."""
    model0.eval()
    model1.eval()
    total_w2_sq = 0.0
    for batch in data:
        Y_pred = make_predicted_barycenter(model0, model1, batch["X0"], batch["X1"], alpha_value)
        total_w2_sq += exact_w2_sq_equal_weights(Y_pred, batch["Y_clean"])
    mean_w2_sq = total_w2_sq / max(len(data), 1)
    rmse_w2 = math.sqrt(max(mean_w2_sq, 0.0))
    return mean_w2_sq, rmse_w2


def evaluate_l2_error(icnn_model, gt_map, bounds=BOUNDS, grid_size=60):
    """Return root mean squared Euclidean map error on a uniform square grid."""
    icnn_model.eval()
    lo, hi = bounds
    grid_val = np.linspace(lo, hi, grid_size, dtype=np.float32)
    xx, yy = np.meshgrid(grid_val, grid_val)
    X_grid = torch.tensor(np.c_[xx.ravel(), yy.ravel()], dtype=torch.float32, device=device)
    Y_gt = gt_map.push_forward_constrained(X_grid).detach().cpu().numpy()
    Y_pred = push_forward_icnn_constrained(icnn_model, X_grid, is_training=False).detach().cpu().numpy()
    mse = np.mean(np.linalg.norm(Y_pred - Y_gt, axis=1) ** 2)
    return float(np.sqrt(mse))


# ------------------------------------------------------------------------------
# Step 5: Validation Stopping and Model-State Management
# ------------------------------------------------------------------------------
class EarlyStopper:
    """Track improvement scaled by max(1, abs(best)) with warm-up and patience."""

    def __init__(self, patience=5, min_rel_improve=1e-3, warmup_epochs=5):
        self.patience = int(patience)
        self.min_rel_improve = float(min_rel_improve)
        self.warmup_epochs = int(warmup_epochs)
        self.best = float("inf")
        self.counter = 0
        self.best_epoch = -1

    def step(self, value, epoch):
        improved = False
        if value < self.best:
            rel_improve = (self.best - value) / max(1.0, abs(self.best))
            if self.best == float("inf") or rel_improve > self.min_rel_improve:
                improved = True

        if improved or epoch < self.warmup_epochs:
            if value < self.best:
                self.best = float(value)
                self.best_epoch = int(epoch)
            self.counter = 0
            return False

        self.counter += 1
        return self.counter >= self.patience


def state_to_cpu(model):
    return {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}


def load_state(model, state):
    model.load_state_dict({k: v.to(device) for k, v in state.items()})


# ------------------------------------------------------------------------------
# Step 6: Single-Replication Block-Coordinate Training
# ------------------------------------------------------------------------------
# Update T0, T1, and alpha0 separately. Keep the other blocks detached so each
# optimizer receives gradients only for its own parameters.
def train_one_trial(trial_idx, N_train, M_particles, alpha_true):
    """Generate one dataset, fit MTDR, and return prediction/parameter diagnostics."""
    seed = int(trial_idx * args.seed_stride)
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # Keep L tied to [-2,2].
    L = (BOUNDS[1] - BOUNDS[0]) / 2.0
    gt_map0 = GroundTruthMap(mode="bubble", strength=0.25, L=L)
    gt_map1 = GroundTruthMap(mode="wave", strength=-0.2, L=L)

    N_val = max(2, int(round(N_train * args.val_ratio)))
    N_test = max(2, int(round(N_train * args.test_ratio)))
    same_cloud_debug = bool(args.same_cloud_debug)

    train_data, val_data, test_data = generate_dataset_splits(
        N_train=N_train,
        N_val=N_val,
        N_test=N_test,
        M_particles=M_particles,
        gt_map0=gt_map0,
        gt_map1=gt_map1,
        alpha_true=alpha_true,
        noise_std=args.noise_std,
        bounds=BOUNDS,
        same_cloud_debug=same_cloud_debug,
    )

    lo, hi = BOUNDS
    if args.anchor_mode == "three":
        anchor_points = [[lo, lo], [hi, lo], [hi, hi]]
    else:
        anchor_points = [[lo, lo], [hi, hi], [lo,hi],[hi,lo]]
    anchors_src = torch.tensor(anchor_points, dtype=torch.float32, device=device)
    anchors_tgt_t0 = gt_map0.push_forward_constrained(anchors_src).detach()
    anchors_tgt_t1 = gt_map1.push_forward_constrained(anchors_src).detach()

    model0 = ICNN(dim=2, hidden_dim=args.hidden_dim).to(device)
    model1 = ICNN(dim=2, hidden_dim=args.hidden_dim).to(device)
    alpha_logit = nn.Parameter(torch.tensor(0.0, dtype=torch.float32, device=device))

    opt0 = optim.Adam(model0.parameters(), lr=args.lr_map)
    opt1 = optim.Adam(model1.parameters(), lr=args.lr_map)
    opt_alpha = optim.Adam([alpha_logit], lr=args.lr_alpha)
    train_loss_fn = SamplesLoss(loss="sinkhorn", p=2, blur=args.blur_train)

    early = EarlyStopper(
        patience=args.early_stop_patience,
        min_rel_improve=args.early_stop_tol,
        warmup_epochs=args.warmup_epochs,
    )

    best_model0_state = state_to_cpu(model0)
    best_model1_state = state_to_cpu(model1)
    best_alpha_logit = float(alpha_logit.detach().cpu().item())
    best_val_w2sq = float("inf")
    best_val_rmse = float("inf")
    final_train_loss = float("nan")
    epochs_run = 0

    for epoch in range(args.outer_epochs):
        epochs_run = epoch + 1
        epoch_loss_sum = 0.0
        epoch_loss_count = 0

        # -----------------------------
        # Phase 1: update T0 only; alpha fixed/detached; T1 fixed.
        # -----------------------------
        model0.train()
        model1.eval()
        for _ in range(args.inner_epochs):
            np.random.shuffle(train_data)
            for batch in train_data:
                opt0.zero_grad(set_to_none=True)
                X0, X1, Y_obs = batch["X0"], batch["X1"], batch["Y_obs"]

                X1_hat_fixed = push_forward_icnn_constrained(model1, X1, is_training=False).detach()
                X0_hat_grad = push_forward_icnn_constrained(model0, X0, is_training=True)

                C = cdist(
                    X0_hat_grad.detach().cpu().numpy(),
                    X1_hat_fixed.cpu().numpy(),
                    metric="sqeuclidean",
                )
                _, col_ind = linear_sum_assignment(C)

                alpha_pred = torch.sigmoid(alpha_logit).detach()
                Y_pred = alpha_pred * X0_hat_grad + (1.0 - alpha_pred) * X1_hat_fixed[col_ind]

                loss_sinkhorn = train_loss_fn(Y_pred, Y_obs)
                T0_anchor_pred = push_forward_icnn_constrained(
                    model0,
                    anchors_src,
                    is_training=True,
                )
                loss_penalty = torch.mean(torch.sum((T0_anchor_pred - anchors_tgt_t0) ** 2, dim=1))
                loss = loss_sinkhorn + args.lambda_penalty * loss_penalty

                loss.backward()
                opt0.step()

                epoch_loss_sum += float(loss.detach().cpu().item())
                epoch_loss_count += 1

        # -----------------------------
        # Phase 2: update T1 only; alpha fixed/detached; T0 fixed.
        # -----------------------------
        model0.eval()
        model1.train()
        for _ in range(args.inner_epochs):
            np.random.shuffle(train_data)
            for batch in train_data:
                opt1.zero_grad(set_to_none=True)
                X0, X1, Y_obs = batch["X0"], batch["X1"], batch["Y_obs"]

                X0_hat_fixed = push_forward_icnn_constrained(model0, X0, is_training=False).detach()
                X1_hat_grad = push_forward_icnn_constrained(model1, X1, is_training=True)

                C = cdist(
                    X1_hat_grad.detach().cpu().numpy(),
                    X0_hat_fixed.cpu().numpy(),
                    metric="sqeuclidean",
                )
                _, col_ind = linear_sum_assignment(C)

                alpha_pred = torch.sigmoid(alpha_logit).detach()
                Y_pred = alpha_pred * X0_hat_fixed[col_ind] + (1.0 - alpha_pred) * X1_hat_grad

                loss_sinkhorn = train_loss_fn(Y_pred, Y_obs)
                T1_anchor_pred = push_forward_icnn_constrained(
                    model1,
                    anchors_src,
                    is_training=True,
                )
                loss_penalty = torch.mean(torch.sum((T1_anchor_pred - anchors_tgt_t1) ** 2, dim=1))
                loss = loss_sinkhorn + args.lambda_penalty * loss_penalty

                loss.backward()
                opt1.step()

                epoch_loss_sum += float(loss.detach().cpu().item())
                epoch_loss_count += 1

        # -----------------------------
        # Phase 3: update alpha only; maps fixed/detached.
        # -----------------------------
        model0.eval()
        model1.eval()
        for _ in range(args.alpha_epochs):
            np.random.shuffle(train_data)
            for batch in train_data:
                opt_alpha.zero_grad(set_to_none=True)
                X0, X1, Y_obs = batch["X0"], batch["X1"], batch["Y_obs"]

                X0_hat = push_forward_icnn_constrained(model0, X0, is_training=False).detach()
                X1_hat = push_forward_icnn_constrained(model1, X1, is_training=False).detach()
                C = cdist(X0_hat.cpu().numpy(), X1_hat.cpu().numpy(), metric="sqeuclidean")
                _, col_ind = linear_sum_assignment(C)

                alpha_pred = torch.sigmoid(alpha_logit)
                Y_pred = alpha_pred * X0_hat + (1.0 - alpha_pred) * X1_hat[col_ind]
                loss_alpha = train_loss_fn(Y_pred, Y_obs)
                loss_alpha.backward()
                opt_alpha.step()

                epoch_loss_sum += float(loss_alpha.detach().cpu().item())
                epoch_loss_count += 1

        final_train_loss = epoch_loss_sum / max(epoch_loss_count, 1)

        # Validation early stopping: exact W2^2, no Sinkhorn blur.
        do_val = ((epoch + 1) % max(args.val_every, 1) == 0) or (epoch == args.outer_epochs - 1)
        if do_val:
            current_alpha = float(torch.sigmoid(alpha_logit).detach().cpu().item())
            val_w2sq, val_rmse = evaluate_prediction_exact(model0, model1, val_data, current_alpha)

            if val_w2sq < best_val_w2sq:
                best_val_w2sq = val_w2sq
                best_val_rmse = val_rmse
                best_model0_state = state_to_cpu(model0)
                best_model1_state = state_to_cpu(model1)
                best_alpha_logit = float(alpha_logit.detach().cpu().item())

            print(
                f"[trial {trial_idx:03d}] epoch {epoch+1:03d}/{args.outer_epochs} | "
                f"alpha={current_alpha:.4f} | train_loss={final_train_loss:.6f} | "
                f"val_W2sq={val_w2sq:.6f} | val_RMSE={val_rmse:.6f} | best_RMSE={best_val_rmse:.6f}",
                flush=True,
            )

            if early.step(val_w2sq, epoch):
                print(
                    f"[trial {trial_idx:03d}] early stop at epoch {epoch+1}; best epoch={early.best_epoch+1}",
                    flush=True,
                )
                break

    # Restore best validation state.
    load_state(model0, best_model0_state)
    load_state(model1, best_model1_state)
    alpha_logit.data = torch.tensor(best_alpha_logit, dtype=torch.float32, device=device)
    final_alpha = float(torch.sigmoid(alpha_logit).detach().cpu().item())

    test_w2sq, test_rmse = evaluate_prediction_exact(model0, model1, test_data, final_alpha)
    l2_err_T0 = evaluate_l2_error(model0, gt_map0, bounds=BOUNDS, grid_size=args.eval_grid)
    l2_err_T1 = evaluate_l2_error(model1, gt_map1, bounds=BOUNDS, grid_size=args.eval_grid)

    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    return {
        "Trial_ID": int(trial_idx),
        "True_Alpha": float(alpha_true),
        "N_Train": int(N_train),
        "N_Val": int(N_val),
        "N_Test": int(N_test),
        "M_Particles": int(M_particles),
        "Est_Alpha": final_alpha,
        "Err_Alpha": abs(final_alpha - alpha_true),
        "Err_T0": l2_err_T0,
        "Err_T1": l2_err_T1,
        "Val_W2sq_Best": best_val_w2sq,
        "Val_RMSE_Best": best_val_rmse,
        "Test_W2sq": test_w2sq,
        "Test_RMSE": test_rmse,
        "Final_Train_Loss": final_train_loss,
        "Epochs_Run": int(epochs_run),
        "Best_Epoch": int(early.best_epoch + 1 if early.best_epoch >= 0 else epochs_run),
        "Bounds_Low": float(BOUNDS[0]),
        "Bounds_High": float(BOUNDS[1]),
        "Blur_Train": float(args.blur_train),
        "LR_Map": float(args.lr_map),
        "LR_Alpha": float(args.lr_alpha),
        "Lambda_Penalty": float(args.lambda_penalty),
        "Anchor_Mode": str(args.anchor_mode),
    }


def worker_task(trial_idx, n_samples, m_particles, alpha_true):
    torch.set_num_threads(1)
    return train_one_trial(trial_idx, n_samples, m_particles, alpha_true)


# ------------------------------------------------------------------------------
# Step 7: Monte Carlo Replications and CSV Output
# ------------------------------------------------------------------------------
# Raw results are checkpointed after each successful replication. Pandas sample
# SDs use ddof=1; with one successful trial the reported SD is NaN.
def main():
    """Run the requested configuration and append its Monte Carlo summary."""
    os.makedirs(args.outdir, exist_ok=True)
    torch.backends.cudnn.benchmark = True
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        pass

    if args.max_workers is None:
        max_workers = 1 if torch.cuda.is_available() else min(8, os.cpu_count() or 1)
    else:
        max_workers = int(args.max_workers)

    if torch.cuda.is_available() and max_workers > 1:
        print(
            "[Warning] You are using CUDA with max_workers > 1. "
            "This can easily cause GPU memory contention. max_workers=1 is recommended on a single GPU.",
            flush=True,
        )

    csv_filename = os.path.join(args.outdir, f"mtdr_icnn_alpha{args.alpha}_N{args.N}_M{args.M}.csv")
    summary_csv = os.path.join(args.outdir, "summary_all.csv")
    results = []
    start_time = time.time()

    print(f">>> Monte Carlo workers: {max_workers}", flush=True)

    def record_result(res):
        results.append(res)
        df_current = pd.DataFrame(results).sort_values("Trial_ID")
        df_current.to_csv(csv_filename, index=False)
        print(
            f"DONE trial {res['Trial_ID']+1:03d}/{args.trials} | "
            f"alpha_hat={res['Est_Alpha']:.4f}, alpha_err={res['Err_Alpha']:.4f}, "
            f"T0={res['Err_T0']:.4f}, T1={res['Err_T1']:.4f}, "
            f"Test_RMSE={res['Test_RMSE']:.6f}, epochs={res['Epochs_Run']}",
            flush=True,
        )

    if max_workers == 1:
        for t in range(args.trials):
            try:
                record_result(worker_task(t, args.N, args.M, args.alpha))
            except Exception as e:
                print(f"ERROR trial {t}: {repr(e)}", flush=True)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_trial = {
                executor.submit(worker_task, t, args.N, args.M, args.alpha): t
                for t in range(args.trials)
            }
            for future in concurrent.futures.as_completed(future_to_trial):
                t = future_to_trial[future]
                try:
                    record_result(future.result())
                except Exception as e:
                    print(f"ERROR trial {t}: {repr(e)}", flush=True)

    if not results:
        raise RuntimeError("No successful trials. Check error messages above.")

    df = pd.DataFrame(results)
    summary = {
        "True_Alpha": args.alpha,
        "N_Train": args.N,
        "M_Particles": args.M,
        "Trials_Success": len(df),
        "Alpha_Est_Mean": df["Est_Alpha"].mean(),
        "Alpha_Est_Std": df["Est_Alpha"].std(),
        "Alpha_Err_Mean": df["Err_Alpha"].mean(),
        "Alpha_Err_Std": df["Err_Alpha"].std(),
        "T0_Err_Mean": df["Err_T0"].mean(),
        "T0_Err_Std": df["Err_T0"].std(),
        "T1_Err_Mean": df["Err_T1"].mean(),
        "T1_Err_Std": df["Err_T1"].std(),
        "Test_W2sq_Mean": df["Test_W2sq"].mean(),
        "Test_W2sq_Std": df["Test_W2sq"].std(),
        "Test_RMSE_Mean": df["Test_RMSE"].mean(),
        "Test_RMSE_Std": df["Test_RMSE"].std(),
        "Val_RMSE_Best_Mean": df["Val_RMSE_Best"].mean(),
        "Final_Train_Loss_Mean": df["Final_Train_Loss"].mean(),
        "Epochs_Run_Mean": df["Epochs_Run"].mean(),
        "Best_Epoch_Mean": df["Best_Epoch"].mean(),
        "Bounds_Low": BOUNDS[0],
        "Bounds_High": BOUNDS[1],
        "Blur_Train": args.blur_train,
        "LR_Map": args.lr_map,
        "LR_Alpha": args.lr_alpha,
        "Lambda_Penalty": args.lambda_penalty,
        "Anchor_Mode": args.anchor_mode,
    }
    summary_df = pd.DataFrame([summary])

    if os.path.exists(summary_csv):
        summary_df.to_csv(summary_csv, mode="a", header=False, index=False)
    else:
        summary_df.to_csv(summary_csv, index=False)

    elapsed_min = (time.time() - start_time) / 60.0
    print(
        f"\n=== Finished {len(df)}/{args.trials} trials in {elapsed_min:.2f} minutes ===",
        flush=True,
    )
    print(summary_df.to_string(index=False), flush=True)
    print(f"\nSaved trial results to: {csv_filename}", flush=True)
    print(f"Saved summary to: {summary_csv}", flush=True)


if __name__ == "__main__":
    main()
