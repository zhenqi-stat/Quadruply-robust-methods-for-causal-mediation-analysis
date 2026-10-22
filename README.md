# Reproducibility code

Code for the numerical studies in *Quadruply robust methods for causal mediation analysis*.

## Experiments

| Script | Experiment | Configuration |
| --- | --- | --- |
| `simu_1.R` | Main text, Setting (i) | `N = 600, 1000, 1400, 1800` |
| `simu_2.R` | Main text, Setting (ii) | `N = 600, 1000, 1400, 1800` |
| `simu_3.R` | Main text, Setting (iii): correctly specified parametric models | `N = 600, 1000, 1400, 1800` |
| `simu_4.R` | Supplement, Setting (iv): binary mediator | `N = 600, 1000, 1400, 1800` |
| `simu_5.R` | Supplement, Setting (v): decayed signal | `signal = 0.25, 0.5, 1, 2`; `N = 1800` |
| `simu_6.R` | Supplement, Setting (vi): varying sparsity | `s_level = 1, 3, 5, 6`; `N = 1800` |
| `real.R` | ACTG175 real-data analysis | 100 sample splits |

The four sparsity levels correspond to `(s_pi, s_q) = (4,12), (7,18), (13,30), (17,39)`, respectively.

## Installation

Run the following command in R:

```r
source("install_packages.R")
```

The required R packages are listed in `install_packages.R`. The ACTG175 dataset is provided by `speff2trial`.

## Running the analyses

Set the working directory to this folder. Select the configuration near the top of the corresponding script, then run it in R:

```r
source("simu_1.R")
```

or from a terminal:

```bash
Rscript simu_1.R
```

All analysis scripts run repetitions `1:100`. MQR1, MQR2, QR1, and QR2 are enabled. Simulation results report Oracle, TR, O-TR, QR1, QR2, O-TR', MQR1, and MQR2, using the manuscript's method names.

For the real-data analysis, run `source("real.R")`.

The entrypoints save results to a new folder under `output/`, including a named `Final_*.RData` file and a summary CSV. Simulation summaries use:

```r
theta <- mean(tau_all)
bias <- apply(pred_theta - theta, 2, median)
RMSE <- sqrt(apply((pred_theta - theta)^2, 2, mean))
```

AL is the median confidence-interval length; AC is empirical coverage. Real-data results aggregate estimates over 100 sample splits.

## Final results

`results/` contains the 25 final RData files used in the manuscript tables: four configurations for each of the six simulation settings, plus `real.RData`. Simulation filenames match their configuration scripts, for example `simu_1_N1800.RData`, `simu_5_signal0.25.RData`, and `simu_6_5.RData`.

Each simulation file contains the estimates and influence-function variances for all eight reported methods over repetitions `1:100`, together with `table_summary`. `real.RData` contains the 100 sample splits for the five methods shown in the real-data table. Source files, method columns, and repetition mappings are recorded within the result files.

To inspect a saved table in R:

```r
result <- new.env()
load("results/simu_1_N1800.RData", envir = result)
result$table_summary
```

## Computing environment

The environment was recorded from the experiment-server R session on 2026-10-22.

| Component | Configuration |
| --- | --- |
| OS | Ubuntu 20.04.4 LTS |
| Processor | AMD EPYC 9654 |
| R | 4.5.2 |
| ranger | 0.18.0 |
| randomForestSRC | 3.9.0 |
| caret | 7.0-1 |
| glmnet | 5.1 |
| BLAS | OpenBLAS 0.3.8 |

Additional package versions, compiler information, and library linkage details are provided in [environment.txt](environment.txt) and [package_versions.csv](package_versions.csv).

## File organization and verification

`configurations/` contains the configuration-specific implementations, selected using `configuration_index.rds`. `source_manifest.csv` maps them to the corresponding result files. `tools/` contains execution, verification, and environment-recording utilities. The final files in `results/` also supply the numerical reference values used by `verify.R`.

Run the archive checks with:

```bash
Rscript verify.R
Rscript verify.R --dgp
```

The first command checks the code and recomputes the table summaries from all 25 saved result files. The second command also regenerates the 24 DGP configurations (four parameter configurations for each of the six simulation settings) and checks each target and 100 Oracle estimates and variances against the references. These checks validate the data-generating processes and Oracle calculations; they do not rerun the nuisance estimation procedures.
