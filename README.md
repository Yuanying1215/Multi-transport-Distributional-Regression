# Multi-transport Distributional Regression (MTDR)

This repository contains the reproduction code for the paper **"Multi-transport Distributional Regression"**.

The proposed framework aggregates predictor-specific transported distributions through a weighted Fréchet mean in the Wasserstein space. It assigns interpretable weights to multiple distributional predictors and defines a flexible regression operator invariant to auxiliary construction choices.

## 1. Prerequisites

The code is written in **R**. Please ensure you have the following packages installed:

``` r
install.packages(c("parallel", "pracma", "fdapace" , "fdadensity" , "quadprog"))
```

## 2. File Structure

The repository is organized as follows:
```text
├── Real_data/               # Code for the real data application (Section 7)
│   ├── Functions.R          # Helper functions for the real data analysis
│   ├── MortFemale.RData     # Pre-processed mortality data (Female)
│   ├── MortMale.RData       # Pre-processed mortality data (Male)
│   └── Results_figure.R     # Script to reproduce the figures and Table 5
│
├── Simu/                    # Code for Simulation Studies (Section 6)
│   ├── 6.1/                 # Single Predictor Simulations (Section 6.1)
│   │   ├── Table1.R         # Reproduces Table 1 (Parameter Estimation)
│   │   └── Table2.R         # Reproduces Table 2 (Comparison: MTDR vs OT vs GOT)
│   │
│   ├── 6.2/                 # Multiple Predictors Simulations (Section 6.2)
│   │   ├── Table3.R         # Reproduces Table 3 (Parameter Estimation, p=2)
│   │   └── Table4.R         # Reproduces Table 4 (Comparison with GOT, p=2)
│   │
│   └── Supp/                # Supplementary Simulations
│       ├── TableS1.R        # Supplementary results for GOT setting
│       └── TableS2.R        # Supplementary results for multi-predictor GOT setting
│
└── README.md
```

## 3. Usage & Reproduction

### Simulation Studies (Section 6)

To reproduce the simulation tables reported in the paper, run the corresponding scripts in the `Simu/` folder.

-   Table 1 (Single Predictor Estimation):

``` bash
Rscript Simu/6.1/Table1.R
```

-   Table 2 (Single Predictor Comparison):

``` bash
Rscript Simu/6.1/Table2.R
```

-   Table 3 (Multiple Predictors Estimation):

``` bash
Rscript Simu/6.2/Table3.R
```

-   Table 4 (Multiple Predictors Comparison):

``` bash
Rscript Simu/6.2/Table4.R
```

### Real Data Application (Section 7)

To reproduce the mortality forecasting results (Table 5) and plots:

``` bash
Rscript Real_data/Results_figure.R
```

## 4. Note on Parallel Computing
The simulation scripts utilize the `parallel` package to accelerate Monte Carlo replications.

* Configuration: Please adjust the number of cores in `makeCluster(20)` within the scripts according to your machine's specifications.

* Default: By default, the scripts utilize `FORK` clusters (supported on Linux/macOS).

* Windows Users: Please modify the `makeCluster` command in the scripts to use `PSOCK` or run sequentially, as `FORK` is not supported on Windows.

