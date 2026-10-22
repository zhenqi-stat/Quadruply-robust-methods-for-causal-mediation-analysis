# Sourcing this file defines the validator only. It never runs fitted methods.
verify_archive_dgp <- function(archive_root) {
  archive_root <- normalizePath(archive_root, mustWork = TRUE)
  index_file <- file.path(archive_root, "configuration_index.rds")
  stopifnot(file.exists(index_file))
  cases <- readRDS(index_file)$cases
  # The bundled final results already contain the target and Oracle references.
  references <- lapply(cases, function(case) {
    path <- file.path(archive_root, "results", paste0(case$id, ".RData"))
    stopifnot(file.exists(path))
    result <- new.env(parent = baseenv())
    load(path, envir = result)
    stopifnot(identical(result$replication_ids, 1:100),
              identical(dim(result$pred_theta), c(100L, 8L)),
              identical(dim(result$pred_sig2), c(100L, 8L)),
              "Oracle" %in% colnames(result$pred_theta),
              "Oracle" %in% colnames(result$pred_sig2))
    list(N = result$N, seed = result$seed, tmax = result$tmax,
         theta = result$theta,
         oracle_theta = unname(result$pred_theta[, "Oracle"]),
         oracle_sig2 = unname(result$pred_sig2[, "Oracle"]))
  })
  stopifnot(length(cases) == 24L, setequal(names(cases), names(references)))

  packages <- c("truncnorm", "mvtnorm", "MASS", "sn")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("DGP verification requires these installed R packages: ",
         paste(missing_packages, collapse = ", "), ". No estimation packages are needed.")
  }

  # Preserve the caller's RNG state; each DGP keeps its own original seed.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  skipped_heads <- c(
    "rm", "library", "require", "source", "sys.source",
    "suppressPackageStartupMessages", "Sys.setenv", "Sys.unsetenv",
    "setwd", "options"
  )

  run_one <- function(case, reference, key) {
    configuration_file <- file.path(archive_root, case$path)
    stopifnot(file.exists(configuration_file))
    lines <- readLines(configuration_file, warn = FALSE)
    expressions <- parse(text = lines, keep.source = TRUE)
    env <- new.env(parent = asNamespace("stats"))
    env$rtruncnorm <- truncnorm::rtruncnorm
    env$rmvnorm <- mvtnorm::rmvnorm
    env$dmvnorm <- mvtnorm::dmvnorm
    env$mvrnorm <- MASS::mvrnorm
    env$rsn <- sn::rsn
    evaluated <- 0L
    last_line <- NA_integer_
    for (i in seq_along(expressions)) {
      expression <- expressions[[i]]
      head <- if (is.call(expression) && is.symbol(expression[[1L]])) as.character(expression[[1L]]) else ""
      if (head %in% skipped_heads) next
      eval(expression, envir = env)
      evaluated <- evaluated + 1L
      # Stop before tuning grids, trainControl, workers, or any nuisance fits.
      if (exists("Y_all", env, inherits = FALSE) && exists("theta", env, inherits = FALSE)) {
        source_reference <- attr(expressions, "srcref")[[i]]
        if (!is.null(source_reference)) last_line <- as.integer(source_reference[[3L]])
        break
      }
    }
    stopifnot(exists("Y_all", env, inherits = FALSE),
              exists("theta", env, inherits = FALSE),
              exists("tau_all", env, inherits = FALSE),
              env$N == case$N, env$N == reference$N,
              env$seed == case$seed, env$seed == reference$seed,
              env$tmax == 100, case$tmax == 100, reference$tmax == 100,
              abs(case$theta - reference$theta) <= 1e-12,
              length(reference$oracle_theta) == 100L,
              length(reference$oracle_sig2) == 100L)

    get_vector <- function(name) as.numeric(get(name, envir = env, inherits = FALSE))
    f0_name <- if (exists("f0_all", env, inherits = FALSE)) "f0_all" else "f_0_all"
    f1_name <- if (exists("f1_all", env, inherits = FALSE)) "f1_all" else "f_1_all"
    A <- get_vector("A_all")
    Y <- get_vector("Y_all")
    pi <- get_vector("pi_all")
    mu <- get_vector("mu_all")
    tau <- get_vector("tau_all")
    f0 <- get_vector(f0_name)
    f1 <- get_vector(f1_name)
    stopifnot(all(lengths(list(A, Y, pi, mu, tau, f0, f1)) == env$N * 100L))

    # Literal row-10 Oracle score from the configuration files. Each column
    # contains one replication's N observations; variance uses divisor N.
    score <- A * f0 / f1 / pi * (Y - mu) +
      (1 - A) / (1 - pi) * (mu - tau) + tau
    score <- matrix(score, nrow = env$N, ncol = 100L)
    oracle_theta <- colMeans(score)
    oracle_sig2 <- colMeans((score - rep(oracle_theta, each = env$N))^2)
    expected_theta <- reference$oracle_theta
    expected_sig2 <- reference$oracle_sig2
    stopifnot(all(is.finite(c(oracle_theta, oracle_sig2, expected_theta, expected_sig2))),
              all(expected_sig2 > 0))
    target_error <- abs(mean(tau) - reference$theta)
    theta_error <- max(abs(oracle_theta - expected_theta))
    sig2_error <- max(abs(oracle_sig2 - expected_sig2))
    sig2_relative_error <- max(abs(oracle_sig2 - expected_sig2) / abs(expected_sig2))
    pass <- target_error <= 1e-12 && theta_error <= 1e-10 &&
      sig2_error <= 1e-8 && sig2_relative_error <= 1e-10
    data.frame(
      case = key, group = case$group, config = case$config,
      N = env$N, seed = env$seed, tmax = env$tmax,
      mean_tau_all = mean(tau), theta_saved = reference$theta,
      theta_abs_difference = target_error,
      oracle_pred_theta_max_abs_difference = theta_error,
      oracle_pred_sig2_max_abs_difference = sig2_error,
      oracle_pred_sig2_max_relative_difference = sig2_relative_error,
      oracle_replications = 100L,
      f0_object = f0_name, f1_object = f1_name,
      last_evaluated_source_line = last_line,
      evaluated_top_level_expressions = evaluated,
      pass = pass, stringsAsFactors = FALSE
    )
  }

  results <- vector("list", length(cases))
  names(results) <- names(cases)
  for (key in names(cases)) {
    message("DGP and Oracle verification: ", key)
    results[[key]] <- tryCatch(
      run_one(cases[[key]], references[[key]], key),
      error = function(error) stop(key, ": ", conditionMessage(error), call. = FALSE)
    )
    gc(verbose = FALSE)
  }
  checks <- do.call(rbind, results)
  rownames(checks) <- NULL
  stopifnot(nrow(checks) == 24L,
            all(checks$theta_abs_difference <= 1e-12),
            all(checks$oracle_pred_theta_max_abs_difference <= 1e-10),
            all(checks$oracle_pred_sig2_max_abs_difference <= 1e-8),
            all(checks$oracle_pred_sig2_max_relative_difference <= 1e-10),
            all(checks$pass))
  checks
}
