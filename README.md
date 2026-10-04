# Multi-transport Distributional Regression (MTDR)

This repository contains the simulation and real-data analysis code for **"Multi-transport Distributional Regression"**.

MTDR aggregates predictor-specific transported distributions through a weighted Wasserstein Fréchet mean. The repository includes two-dimensional ICNN-based simulations, one-dimensional supplementary simulations and comparisons with OT and GOT, and a mortality-data application with reference-distribution sensitivity analysis.

## 1. Prerequisites

The main simulation scripts use **Python**. The supplementary simulations and real-data application use **R**.

### Python

Install the dependencies in your Python environment:

```bash
python -m pip install numpy pandas scipy torch geomloss
```

The scripts use CUDA when available and otherwise run on the CPU. GPU execution requires a PyTorch installation compatible with your CUDA environment. Numerical settings and command-line options are documented in each script.

### R

```r
install.packages(c("pracma", "fdapace", "fdadensity", "frechet"))
```

The `parallel` package is included with R. The `frechet` package is used by the density plots in `Results_figure.R`; reference sensitivity itself requires only `pracma` and `fdadensity`. The plotting code calls the internal function `frechet:::qf2pdf`, so compatibility may depend on the installed package version.

## 2. File Structure

The repository uses the directory names `Simu` and `Real_data`.

```text
├── Real_data/                  # Mortality-data application
│   ├── Functions.R             # MTDR, OT, and GOT fitting functions
│   ├── MortFemale.RData        # Preprocessed female mortality data
│   ├── MortMale.RData          # Preprocessed male mortality data
│   ├── Results_figure.R        # Model comparisons and Bulgaria plots
│   └── reference_sensitivity.R # MTDR refits under three reference choices
│
├── Simu/
│   ├── 6.1/
│   │   └── Table1.py           # Two-dimensional MTDR, one random predictor
│   ├── 6.2/
│   │   └── Table2.py           # Two-dimensional MTDR, two random predictors
│   └── Supp/                  # One-dimensional supplementary simulations
│       ├── TableS1.R           # Single predictor: estimation and prediction
│       ├── TableS2.R           # MTDR vs. OT vs. GOT under the MTDR setting
│       ├── TableS3.R           # MTDR vs. OT vs. GOT under the GOT setting
│       ├── TableS4.R           # Multiple predictors: estimation and prediction
│       ├── TableS5.R           # MTDR vs. GOT under the MTDR setting
│       └── TableS6.R           # MTDR vs. GOT under the GOT setting
│
└── README.md
```

## 3. Simulation Studies

Run the commands in this section from the **repository root**. The Python commands below each run one configuration, not an entire table. Repeat with the weight and sample-size combinations of interest, using distinct raw-output paths for separate runs.

### Table 1: Two-dimensional distributions, one random predictor

```bash
python Simu/6.1/Table1.py \
    --alpha 0.5 --N 50 --M 400 --trials 50 \
    --gpu_id 0 --max_workers 1 \
    --outdir results/table1_alpha05_N50_M400
```

- `--alpha` is the true **reference weight** `alpha0`; the predictor weight is `1 - alpha0`.
- `--N` is the number of training distributions and `--M` is the number of particles per distribution.
- Validation and test sizes default to `max(2, round(0.2*N))` and `max(2, round(0.3*N))`.
- Each trial uses seed `trial_index * seed_stride`, with default stride 1000.

The output directory is created automatically. It contains a trial-level CSV named `mtdr_icnn_alpha{alpha}_N{N}_M{M}.csv` and an appended `summary_all.csv`. Prediction summaries include `Test_W2sq_Mean` and `Test_RMSE_Mean`, together with their Monte Carlo standard deviations. Evaluation uses exact balanced empirical squared Wasserstein distances via optimal matching; the trial-level RMSE is the square root of the mean squared test distance. Training instead uses a Sinkhorn loss.

### Table 2: Two-dimensional distributions, two random predictors

```bash
mkdir -p results/table2
python Simu/6.2/Table2.py \
    --alpha_true 0.3,0.35,0.35 \
    --N 50 --N_val 15 --N_test 15 --M 200 --trials 50 \
    --gpu_id 0 \
    --out_csv results/table2/alpha030_035_035_N50_M200.csv \
    --summary_csv results/table2/summary_all.csv
```

- `--alpha_true` lists the reference weight followed by the two predictor weights. Nonnegative inputs with positive sum are normalized to sum to one.
- `--N_val` and `--N_test` are explicit counts, not ratios of `--N`.
- The default number of trials is **1**; specify `--trials 50` for 50 replications.
- Replication `t` uses seed `seed + 1000*t`; the base seed can be set with `--seed`.

This script uses fixed-support entropic barycenter computation. The prediction summary is `test_sinkhorn_loss_mean`, with `test_sinkhorn_loss_std` across replications. This criterion is a **Sinkhorn divergence**, without a square root, not the exact matching-based Wasserstein criterion in `Table1.py`. Barycenter regularization (`--epsilon_bary`) and prediction-loss smoothing (`--blur_loss`) are separate settings.

The raw CSV is updated after each successful trial; the summary CSV accumulates rows across runs. Parent output directories must already exist. Neither Python script assembles a publication-ready table automatically.

To inspect all Python options without starting simulation:

```bash
python Simu/6.1/Table1.py --help
python Simu/6.2/Table2.py --help
```

### Tables S1-S6: One-dimensional supplementary simulations

Each R script contains its own data generation, fitting routines, and simulation settings. Run the desired script:

```bash
Rscript Simu/Supp/TableS1.R
Rscript Simu/Supp/TableS2.R
Rscript Simu/Supp/TableS3.R
Rscript Simu/Supp/TableS4.R
Rscript Simu/Supp/TableS5.R
Rscript Simu/Supp/TableS6.R
```

These scripts report Monte Carlo summaries to the console. Sample sizes, weights, replication counts, and worker counts are set inside the scripts rather than through command-line arguments. Check the parallel-computing notes below before running them.

## 4. Real Data Application

The application predicts male age-at-death distributions in 2010 from male and female distributions in 2005 for 34 countries. Keep both `.RData` files and `Functions.R` in `Real_data`.

### Model comparisons and figures

From the repository root, first change to the data directory so that the script can find its relative input paths:

```bash
cd Real_data
Rscript Results_figure.R
```

The script fits MTDR, OT, and GOT, evaluates leave-one-out prediction errors, and plots the Bulgaria example. Alternatively, set the R/RStudio working directory to `Real_data` and run `source("Results_figure.R")` to inspect the results and plots interactively. The script does not export a formatted manuscript table.

### Reference-distribution sensitivity (supplementary Tables S7-S8)

Run from `Real_data`:

```bash
Rscript --vanilla reference_sensitivity.R
```

This script refits **MTDR only**; it does not fit OT or GOT. It compares three references on the age interval `[0, 100]`:

- **FM** (`response_mean`): the Wasserstein Fréchet mean of the training responses in each fold.
- **U** (`uniform`): the uniform distribution on `[0, 100]`.
- **TN** (`truncated_normal`): a normal distribution with mean 50 and standard deviation 25, truncated to `[0, 100]`.

By default, all leave-one-out folds and three additional full-sample fits are run. Results are saved in a new `reference_results_YYYYMMDD_HHMMSS` directory. The principal outputs are:

| File | Contents |
| --- | --- |
| `loo_summary.csv` | Reference-level LOOCV summaries. `prediction_W2_mean` is the average Wasserstein distance (AWD). |
| `loo_by_country.csv` | Held-out prediction errors, fitted weights, cross-reference differences, and fitting diagnostics for each country/reference. |
| `full_sample_weights.csv` | Weights from fitting all countries together, separate from LOOCV; omitted when `--full-fit=false`. |
| `fits/`, `logs/` | Fit checkpoints and optimization logs. |
| `config.rds`, `config.txt`, `sessionInfo.txt`, `prepared_data.rds`, `Functions_used.R` | Settings, software information, prepared data, and a copy of the fitting functions. |

Weight order is `alpha0_reference`, `alpha1_male2005`, `alpha2_female2005`. Weights in a held-out country's row are estimated from the other countries; they are not country-specific regression coefficients. Wasserstein distances here are measured in years.

All sensitivity options use `--name=value`. For example, to skip full-sample fits and choose an output directory:

```bash
Rscript --vanilla reference_sensitivity.R \
    --full-fit=false --output-dir=reference_results_loo
```

To run from the repository root instead, specify the data directory explicitly:

```bash
Rscript --vanilla Real_data/reference_sensitivity.R --data-dir=Real_data
```

The default iteration budget is 500 and the stopping tolerance is `1e-7`. Inspect fit status, warnings, and iteration-limit diagnostics before interpreting the comparisons. To resume an interrupted analysis, supply the same output directory and analysis options with `--resume=true`; saved fits, including flagged fits, are reused rather than rerun.

## 5. Computing and Reproducibility Notes

- **Python simulations:** `Table1.py` supports parallel Monte Carlo workers; `--max_workers 1` is recommended on a single GPU. `Table2.py` runs trials sequentially. Both scripts use validation-based stopping and evaluate on separate test data.
- **R simulations:** the supplementary scripts use `makeCluster(20, type="FORK")`. Adjust 20 to suit the available resources. `FORK` is supported on Linux/macOS, not Windows. A Windows port requires configuring PSOCK workers, including their packages, functions, and data; changing the cluster type alone may not suffice.
- **Outputs:** Python summaries include successful trials only. Check the successful-trial count and logs, and use separate raw-output paths to avoid overwriting earlier runs. Monte Carlo standard deviations require at least two successful replications.
- **Numerical reproducibility:** seeds and settings are recorded or defined in the scripts, but identical results across hardware and package versions are not guaranteed. Retain the commands, output files, and software versions used for each reported experiment.
