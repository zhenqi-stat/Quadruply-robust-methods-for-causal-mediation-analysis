# Data dictionary

## Source data and preparation

The real-data analysis uses `ACTG175` from `speff2trial` 1.0.5. Install the package with `install_packages.R`, then access the original 2,139-row, 27-column data frame in R:

```r
data("ACTG175", package = "speff2trial")
help("ACTG175", package = "speff2trial")
```

Package source and documentation: https://CRAN.R-project.org/package=speff2trial and https://cran.r-project.org/web/packages/speff2trial/speff2trial.pdf.

The code keeps `offtrt == 0`, `r == 1`, and `arms != 3`, preserving the original row order. This yields 780 participants (253 with `treat == 0`, 527 with `treat == 1`). None of the selected analysis variables is missing. No imputation is applied. The source data remain available through the package and are not duplicated in `results/`.

| Source variable | Type and unit or coding | Role |
| --- | --- | --- |
| `age` | Numeric, baseline years | Covariate |
| `wtkg` | Numeric, baseline kg | Covariate |
| `karnof` | Numeric, score 0-100 | Covariate |
| `cd40`, `cd80` | Numeric, baseline CD4/CD8 counts per cubic millimeter | Covariates |
| `hemo`, `homo`, `drugs` | Binary indicators for hemophilia, homosexual activity, IV drug-use history | Covariates |
| `race` | 0 white, 1 non-white | Covariate |
| `gender` | 0 female, 1 male | Covariate |
| `str2` | 0 treatment-naive, 1 experienced | Covariate |
| `symptom` | 0 asymptomatic, 1 symptomatic | Covariate |
| `treat` | 0 zidovudine alone, 1 other treatment | Treatment `A_all`; arm 3 is excluded |
| `arms` | 0 zidovudine, 1 zidovudine/didanosine, 2 zidovudine/zalcitabine, 3 didanosine | Selection |
| `offtrt` | 1 off treatment before 96 +/- 5 weeks, 0 otherwise | Selection |
| `r` | 1 week-96 CD4 observed, 0 missing | Selection |
| `cd420`, `cd820` | Numeric, CD4/CD8 counts at 20 +/- 5 weeks | Mediators, natural-log transformed |
| `cd496` | Numeric, CD4 count at 96 +/- 5 weeks | Outcome, natural-log transformed |

`X_all` contains the 12 listed baseline covariates. `M_all = log(df[, c("cd420", "cd820")])`, `Y_all = log(df$cd496)`, and `A_all = df$treat`. Polynomial and interaction bases are constructed by `expand_covariates_with_mediation` in `configurations/real_2_data_generate.R`. The QR-specific mediator basis is constructed after the log transformation. All preprocessing is executed by `real.R` through `configurations/real_2.R`.

## Simulation inputs

Each `configurations/simu_*.R` file defines a DGP, generates a pool of `N * tmax` observations with `tmax = 100`, and defines `X_all` (covariates), `A_all` (binary treatment), `M_all` (mediators), `Y_all` (outcome), and the true nuisance functions. `tau_all` is the true integrated outcome regression evaluated on that pool. The target is the finite-pool average `theta = mean(tau_all)`. Each replication uses its designated `N`-row block. Changing the executed replication subset does not change `tmax` or the generated pool.

`configuration_index.rds` is an R list with a `cases` component containing the 24 configuration records (`id`, `group`, `config`, `path`, `theta`, `N`, `seed`, `tmax`). `source_manifest.csv` maps configuration implementations to the saved results. `table_map.csv` maps the public entrypoints to the manuscript tables.

## Final simulation results

The 24 simulation files in `results/` are final computational results, not complete raw generated datasets. They contain 100 replications and eight reported methods.

| Object | Meaning |
| --- | --- |
| `source_id`, `source_path`, `group`, `config` | Configuration identifier, relative implementation path, family, selector |
| `N`, `seed`, `tmax`, `theta` | Sample size, DGP seed, pool replication count, target |
| `replication_ids` | Integers 1 through 100; row identifiers |
| `method_names` | Column labels: MQR1, MQR2, O-TR', QR1, QR2, TR, O-TR, Oracle |
| `pred_theta` | 100 x 8 matrix of replication-level estimates |
| `pred_sig2` | 100 x 8 matrix of influence-function variance estimates; variance of the estimator is `pred_sig2 / N` |
| `results_par` | The eight estimate columns followed by the eight variance columns, explicitly named |
| `bias` | Columnwise `median(pred_theta - theta)` |
| `RMSE` | Columnwise `sqrt(mean((pred_theta - theta)^2))` |
| `AL`, `AC` | Median two-sided 95% normal CI length; empirical coverage proportion |
| `ESD`, `ASD` | 1.4826 times the median absolute deviation of estimates; median estimated standard error |
| `sig2` | Median replication-level influence-function variance |
| `table_summary`, `manuscript_table` | Unrounded numerical summaries and rounded manuscript values |
| `source_map` | Replication-by-method provenance: input basename/hash, row, and method column |
| `schema_version`, `theta_definition`, `paper_method_names` | Format version, target definition, manuscript ordering |

Bias, RMSE, and AL are displayed to three decimal places; AC is displayed as a proportion to two decimal places. No multiplier is applied.

## Final real-data results

`results/real.RData` contains five reported methods over 100 sample splits. The columns are MQR2, O-TR', QR2, TR, and O-TR. Matrices `pred_theta`, `pred_theta0`, and `pred_theta1` hold estimates of theta(1,0), theta(0,0), and theta(1,1), respectively. `pred_sig2` and `pred_sig2_DE`, `pred_sig2_IDE`, `pred_sig2_TE` hold the corresponding influence-function variance estimates. NDE, NIE, and ATE splitwise estimates are `pred_theta - pred_theta0`, `pred_theta1 - pred_theta`, and `pred_theta1 - pred_theta0`.

Each point estimate is the median over splits. For each scalar effect, the aggregated influence-function variance is `median(split_variance + (split_estimate - median(split_estimate))^2)`. The standard error is the square root of that quantity divided by `N`. `result_summary` and `table_summary` hold the aggregate estimates, CI lengths, p-values, and `NIE / ATE`; `confidence_intervals` holds the CI endpoints. `manuscript_display` stores the exact table strings. `aggregation` records these definitions, and `source_recipe` records source-file, method-column, and split-range provenance.

`results_par` concatenates seven named split-level matrices. `N = 780`; `replication_ids` identifies splits 1 through 100. `DE`, `IDE`, and `TE` are the aggregate NDE, NIE, and ATE; the `p_*`, `sig2_*`, and `theta*` objects contain their p-values, variances, and component means. The mediation proportion is reported as a percentage in the manuscript. These result files contain no patient-level records.
