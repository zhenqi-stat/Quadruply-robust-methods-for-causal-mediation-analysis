# Collect reproducibility metadata only: no model fits or package installation.
collect_environment <- function(output_parent = ".") {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!dir.exists(output_parent)) dir.create(output_parent, recursive = TRUE)
  output_parent <- normalizePath(output_parent, mustWork = TRUE)
  output_directory <- tempfile(
    pattern = paste0("environment_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
    tmpdir = output_parent
  )
  stopifnot(dir.create(output_directory))

  packages <- c(
    "butcher", "caret", "cobalt", "corpcor", "densratio", "doParallel",
    "doRNG", "dplyr", "filelock", "FNN", "foreach", "glmnet", "itertools",
    "ks", "logger", "MASS", "matrixStats", "mice", "mvtnorm", "purrr",
    "randomForestSRC", "ranger", "RCAL", "sn", "speff2trial", "survival",
    "tidyr", "truncnorm", "Matrix", "Rcpp", "RcppEigen", "pROC"
  )
  installed <- utils::installed.packages()
  package_information <- lapply(packages, function(package) {
    rows <- which(installed[, "Package"] == package)
    if (!length(rows)) return(data.frame(package = package, version = NA_character_,
      built = NA_character_, library = NA_character_))
    data.frame(package = package, version = installed[rows, "Version"],
      built = installed[rows, "Built"], library = installed[rows, "LibPath"],
      row.names = NULL)
  })
  package_information <- do.call(rbind, package_information)
  run_command <- function(command, arguments = character()) {
    path <- unname(Sys.which(command))
    if (!nzchar(path) && !file.exists(command)) return("Command not available.")
    if (!nzchar(path)) path <- command
    tryCatch(suppressWarnings(system2(path, arguments, stdout = TRUE, stderr = TRUE)),
      error = function(error) conditionMessage(error))
  }
  r_binary <- file.path(R.home("bin"), "R")
  r_config_names <- c("CC", "CXX", "CXX17", "CXXFLAGS", "CXX17FLAGS",
                      "LDFLAGS", "BLAS_LIBS", "LAPACK_LIBS", "FLIBS")
  r_config <- setNames(lapply(r_config_names, function(name) {
    run_command(r_binary, c("CMD", "config", name))
  }), r_config_names)
  compilers <- setNames(lapply(c("gcc", "g++", "gfortran", "clang"), function(name) {
    head(run_command(name, "--version"), 3L)
  }), c("gcc", "g++", "gfortran", "clang"))

  ranger_library <- ranger_linkage <- ranger_md5 <- NA_character_
  if (requireNamespace("ranger", quietly = TRUE)) {
    dll <- getLoadedDLLs()[["ranger"]]
    if (!is.null(dll)) {
      ranger_library <- dll[["path"]]
      ranger_md5 <- unname(tools::md5sum(ranger_library))
      if (identical(Sys.info()[["sysname"]], "Linux")) {
        ranger_linkage <- run_command("ldd", shQuote(ranger_library))
      } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
        ranger_linkage <- run_command("otool", c("-L", shQuote(ranger_library)))
      }
    }
  }
  thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "R_RANGER_NUM_THREADS", "RF_CORES")
  record <- list(
    collected_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    system = Sys.info(), R_version = R.version, session_info = utils::sessionInfo(),
    external_software = extSoftVersion(), library_paths = .libPaths(),
    packages = package_information, RNG_kind = RNGkind(),
    logical_cores = parallel::detectCores(logical = TRUE),
    physical_cores = parallel::detectCores(logical = FALSE),
    thread_environment = Sys.getenv(thread_variables, unset = NA_character_),
    thread_options = options()[intersect(names(options()), c("ranger.num.threads", "rf.cores", "mc.cores"))],
    R_build_configuration = r_config, installed_compiler_commands = compilers,
    ranger_library = ranger_library, ranger_library_md5 = ranger_md5,
    ranger_linkage = ranger_linkage,
    note = paste("Installed compiler commands and R configuration do not necessarily identify",
                 "the compiler or flags originally used to build an installed ranger binary.")
  )
  if (identical(Sys.info()[["sysname"]], "Linux")) {
    record$cpu <- run_command("lscpu")
    record$os_release <- if (file.exists("/etc/os-release")) readLines("/etc/os-release", warn = FALSE) else NA_character_
    record$runtime_package_builds <- run_command("dpkg-query",
      c("-W", "libstdc++6", "libc6", "libgomp1", "libopenblas0-pthread"))
  } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
    record$cpu <- run_command("sysctl", c("-n", "machdep.cpu.brand_string"))
  }
  utils::write.csv(package_information, file.path(output_directory, "package_versions.csv"), row.names = FALSE)
  saveRDS(record, file.path(output_directory, "environment.rds"))
  capture.output(print(record), file = file.path(output_directory, "environment.txt"))
  message("Environment metadata saved to: ", output_directory)
  invisible(output_directory)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  collect_environment(if (length(arguments)) arguments[[1L]] else ".")
}
