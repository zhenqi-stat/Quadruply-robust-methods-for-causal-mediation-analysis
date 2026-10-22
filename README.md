# Reproducibility code

Code and final computational results for *Quadruply robust methods for causal mediation analysis*.

## Experiments and numerical tables

| Script | Experiment | Table | Configurations |
| --- | --- | --- | --- |
| `simu_1.R` | Setting (i) | Main Table 5 | `N = 600, 1000, 1400, 1800` |
| `simu_2.R` | Setting (ii) | Main Table 6 | `N = 600, 1000, 1400, 1800` |
| `simu_3.R` | Setting (iii), correctly specified parametric models | Main Table 7 | `N = 600, 1000, 1400, 1800` |
| `real.R` | ACTG175 analysis | Main Table 8 | 100 sample splits |
| `simu_4.R` | Setting (iv), binary mediator | Supplement Table 4 | `N = 600, 1000, 1400, 1800` |
| `simu_5.R` | Setting (v), decayed signal | Supplement Table 5 | `signal = 0.25, 0.5, 1, 2`; `N = 1800` |
| `simu_6.R` | Setting (vi), varying sparsity | Supplement Table 6 | `s_level = 1, 3, 5, 6`; `N = 1800` |

The four sparsity levels correspond to `(s_pi, s_q) = (4,12), (7,18), (13,30), (17,39)`. `table_map.csv` gives the same mapping in machine-readable form. The other manuscript tables summarize theoretical properties, and the figures are schematic diagrams.

## Installation and data

Set the R working directory to this folder and run:

```r
source("install_packages.R")
```

The required packages are listed in `install_packages.R`; recorded versions are in `package_versions.csv`. ACTG175 is distributed with `speff2trial` 1.0.5 at https://CRAN.R-project.org/package=speff2trial. It is loaded automatically by the real-data script. `DATA_DICTIONARY.md` documents the analysis variables, participant selection, transformations, and result-file objects.

## Regenerate the numerical tables

```r
source("reproduce_tables.R")
```

or, from a terminal:

```bash
Rscript run_all.R
```

This recalculates the seven numerical tables from the 25 final result files in `results/` and writes unrounded CSVs, formatted CSVs, and LaTeX rows to `tables/`. It takes less than one minute and requires only base R. The source result files are not overwritten.

## Rerun estimation

Select the configuration near the top of an entrypoint, then run it, for example:

```r
source("simu_1.R")
```

Run the full collection of 24 simulation configurations and the real-data analysis with:

```bash
Rscript run_all.R --fit
```

All analysis scripts run repetitions `1:100`. MQR1, MQR2, QR1, and QR2 are enabled. Simulation results report Oracle, TR, O-TR, QR1, QR2, O-TR', MQR1, and MQR2. The real-data table uses TR, O-TR, QR2, O-TR', and MQR2; the real-data implementation also computes QR1 and MQR1.

Simulations use 30 worker processes. The real-data implementation specifies 31 workers. Its QR2, TR, and O-TR batch took about 2 hours for 100 sample splits. Approximate runtimes and batch coverage are summarized in `RUNTIME.md`.

### Progress and output files

The run directory is printed when execution starts. Method-by-method progress is written to log files, not continuously to the R console or terminal. A quiet console during model fitting is normal. Look under the newly created directory:

```text
output/<configuration>/<run>/internal_raw/logs/     Settings (i), (ii), (iii), (v), (vi)
output/<configuration>/<run>/internal_raw/ld1/     Setting (iv)
output/real/<run>/internal_raw/logs_cd/             ACTG175
```

Each completed run writes a named `Final_*.RData`, a summary CSV, and `run_manifest.rds` in its run directory. `internal_raw/` holds the implementation's intermediate outputs and logs. New runs do not replace the bundled final results.

## Summary definitions

Simulation summaries use:

```r
theta <- mean(tau_all)
bias <- apply(pred_theta - theta, 2, median)
RMSE <- sqrt(apply((pred_theta - theta)^2, 2, mean))
```

AL is the median 95% confidence-interval length; AC is empirical coverage. Real-data results aggregate estimates over 100 sample splits. Full definitions are in `DATA_DICTIONARY.md` and `tools/verify_results.R`.

## Computing environment

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

Additional package versions, compiler information, and library-linkage details are provided in `environment.txt` and `package_versions.csv`. `SOFTWARE_REFERENCES.bib` contains references for the principal packages.

## Structure and verification

`configurations/` contains the configuration-specific implementations. `configuration_index.rds` selects their parameters; `source_manifest.csv` and the result-file metadata record source and replication mappings. `tools/` contains execution and verification utilities. `results/` holds 24 simulation result files, each with 100 replications and eight reported methods, plus `real.RData` with 100 splits for the five methods in the real-data table.

The configuration-specific files retain their separate implementations so that their tuning rules and random-number call order stay tied to the corresponding results. Inactive alternative snippets in comments are not part of execution. Internal companion calculations are retained where required to preserve the random-number sequence; their results are not included among the reported methods.

```bash
Rscript verify.R
Rscript verify.R --dgp
```

The first command checks the source files, result completeness, and table calculations. The second also regenerates the 24 DGP configurations and checks each finite-pool target and 100 Oracle estimates and variances. These DGP/Oracle checks do not rerun nuisance estimation.

## License

The supplied code is released under the MIT License; see `LICENSE`. Third-party R packages and the ACTG175 data retain their respective terms.
