# Install missing direct dependencies; do not upgrade an existing environment.
# Run explicitly with source("install_packages.R") or Rscript install_packages.R.
simulation_packages <- c(
  "butcher", "caret", "cobalt", "corpcor", "densratio", "doParallel",
  "filelock", "FNN", "foreach", "glmnet", "itertools", "logger", "MASS",
  "matrixStats", "mvtnorm", "randomForestSRC", "ranger", "RCAL", "sn",
  "truncnorm"
)
real_additional_packages <- c(
  "doRNG", "dplyr", "ks", "mice", "purrr", "speff2trial", "survival", "tidyr"
)
required_packages <- unique(c(simulation_packages, real_additional_packages))
missing_packages <- setdiff(required_packages, rownames(utils::installed.packages()))
if (length(missing_packages)) {
  utils::install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
} else {
  message("All direct dependencies are already installed; no packages were changed.")
}
remaining <- setdiff(required_packages, rownames(utils::installed.packages()))
if (length(remaining)) {
  stop("These dependencies are still missing: ", paste(remaining, collapse = ", "))
}
# CRAN resolves current compatible releases, not the recorded historical builds.
