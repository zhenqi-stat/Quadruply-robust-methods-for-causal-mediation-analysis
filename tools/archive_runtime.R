# Portable execution and clean output for the retained reviewer methods.
# The configuration scripts contain the scientific implementation. This helper
# only selects a script, controls the requested replication IDs, and exports it.

.archive_stop <- function(...) stop(..., call. = FALSE)

.archive_replication_ids <- function(replication_ids) {
  if (!is.numeric(replication_ids) || !length(replication_ids) ||
      anyNA(replication_ids) || any(!is.finite(replication_ids)) ||
      any(replication_ids != floor(replication_ids)) ||
      any(replication_ids < 1 | replication_ids > 100) ||
      anyDuplicated(replication_ids)) {
    .archive_stop("replication_ids must contain distinct integer IDs from 1 through 100.")
  }
  as.integer(replication_ids)
}

.archive_cases <- function(index) {
  cases <- index$cases
  if (is.data.frame(cases)) {
    return(lapply(seq_len(nrow(cases)), function(i) {
      lapply(cases[i, , drop = FALSE], function(x) {
        if (is.list(x)) x[[1L]] else x[[1L]]
      })
    }))
  }
  if (!is.list(cases) || !length(cases)) {
    .archive_stop("configuration_index.rds must contain a nonempty cases list or data frame.")
  }
  cases
}

.archive_select_case <- function(archive_root, group, config) {
  index_path <- file.path(archive_root, "configuration_index.rds")
  if (!file.exists(index_path)) .archive_stop("Missing configuration index: ", index_path)
  cases <- .archive_cases(readRDS(index_path))
  matches <- vapply(cases, function(case) {
    identical(as.character(case$group), as.character(group)) &&
      identical(as.character(case$config), as.character(config))
  }, logical(1))
  if (sum(matches) != 1L) {
    .archive_stop("Expected one configuration for simulation ", group,
                  " with selector ", config, "; found ", sum(matches), ".")
  }
  case <- cases[[which(matches)]]
  required <- c("id", "group", "config", "path", "theta", "N", "seed", "tmax")
  if (!all(required %in% names(case))) {
    .archive_stop("The selected configuration lacks: ",
                  paste(setdiff(required, names(case)), collapse = ", "))
  }
  if (length(case$tmax) != 1L || !isTRUE(case$tmax == 100)) {
    .archive_stop("The scientific configuration must retain tmax = 100.")
  }
  case
}

.archive_source_path <- function(archive_root, path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .archive_stop("Each configuration path must be one relative path.")
  }
  root <- normalizePath(archive_root, mustWork = TRUE)
  candidate <- normalizePath(file.path(root, path), mustWork = TRUE)
  if (!startsWith(candidate, paste0(root, .Platform$file.sep))) {
    .archive_stop("The configuration must be inside the archive: ", path)
  }
  candidate
}

.archive_lhs_is <- function(expression, name) {
  is.call(expression) && length(expression) >= 3L &&
    is.symbol(expression[[1L]]) &&
    as.character(expression[[1L]]) %in% c("<-", "=", "<<-") &&
    is.symbol(expression[[2L]]) && identical(as.character(expression[[2L]]), name)
}

.archive_output_dimensions <- function(expression) {
  if (!is.call(expression) || identical(expression[[1L]], as.name("function"))) {
    return(expression)
  }
  for (i in seq_along(expression)) {
    if (identical(expression[[i]], quote(expr = ))) next
    expression[[i]] <- .archive_output_dimensions(expression[[i]])
  }
  if (length(expression) >= 4L && identical(expression[[1L]], as.name("[")) &&
      identical(expression[[2L]], as.name("results_par")) &&
      !"drop" %in% names(expression)) {
    expression$drop <- FALSE
  }
  expression
}

.archive_evaluate <- function(path, replication_ids, execution, real = FALSE) {
  expressions <- parse(path, keep.source = TRUE)
  found_replication_ids <- FALSE
  collected_results <- FALSE
  for (expression in expressions) {
    if (collected_results) expression <- .archive_output_dimensions(expression)
    if (.archive_lhs_is(expression, "replication_ids")) {
      assign("replication_ids", replication_ids, envir = execution)
      found_replication_ids <- TRUE
    } else if (real && .archive_lhs_is(expression, "v") &&
               identical(expression[[3L]], as.name("tmax"))) {
      # v controls aggregation over completed splits; tmax remains the DGP size.
      assign("v", length(replication_ids), envir = execution)
    } else {
      eval(expression, envir = execution)
    }
    if (.archive_lhs_is(expression, "results_par")) {
      value <- get("results_par", envir = execution, inherits = FALSE)
      if (length(replication_ids) == 1L && is.numeric(value) && is.null(dim(value))) {
        assign("results_par", matrix(value, nrow = 1L), envir = execution)
      }
      collected_results <- TRUE
    }
  }
  if (!found_replication_ids) {
    .archive_stop("The prepared configuration must define replication_ids before its foreach loop.")
  }
  if (!exists("tmax", envir = execution, inherits = FALSE) ||
      !isTRUE(get("tmax", envir = execution) == 100)) {
    .archive_stop("The configuration changed tmax; it must remain 100 for scientific reproducibility.")
  }
  invisible(execution)
}

.archive_matrix <- function(execution, name, n, p) {
  if (!exists(name, envir = execution, inherits = FALSE)) {
    .archive_stop("The configuration did not produce ", name, ".")
  }
  value <- get(name, envir = execution, inherits = FALSE)
  # foreach returns one row as a vector for a one-replication execution.
  if (is.numeric(value) && is.null(dim(value)) && n == 1L && length(value) == p) {
    value <- matrix(value, nrow = 1L)
  }
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(as.integer(dim(value)), as.integer(c(n, p)))) {
    .archive_stop(name, " must be a numeric ", n, " by ", p,
                  " matrix; the replication range must not change tmax.")
  }
  value
}

.archive_require_finite <- function(value, name) {
  if (is.numeric(value) && any(!is.finite(value))) {
    .archive_stop("Non-finite primary output in ", name,
                  ". All retained methods must be enabled and complete before clean export.")
  }
  if (is.data.frame(value) || is.list(value)) {
    for (field in names(value)) .archive_require_finite(value[[field]], paste0(name, "$", field))
  }
  invisible(TRUE)
}

.archive_simulation_result <- function(execution, case, replication_ids) {
  n <- length(replication_ids)
  historical_columns <- c(1L, 2L, 3L, 5L, 6L, 7L, 8L, 10L)
  method_names <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
  pred_theta <- .archive_matrix(execution, "pred_theta", n, 10L)[, historical_columns, drop = FALSE]
  pred_sig2 <- .archive_matrix(execution, "pred_sig2", n, 10L)[, historical_columns, drop = FALSE]
  dimnames(pred_theta) <- dimnames(pred_sig2) <- list(as.character(replication_ids), method_names)
  .archive_require_finite(pred_theta, "pred_theta")
  .archive_require_finite(pred_sig2, "pred_sig2")
  if (any(pred_sig2 < 0)) .archive_stop("Retained influence-function variances must be nonnegative.")
  theta <- get("theta", envir = execution, inherits = FALSE)
  N <- get("N", envir = execution, inherits = FALSE)
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
      !is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    .archive_stop("Simulation theta must be finite and N must be a positive finite scalar.")
  }
  if (!isTRUE(all.equal(as.numeric(theta), as.numeric(case$theta), tolerance = 1e-12)) ||
      !isTRUE(all.equal(as.numeric(N), as.numeric(case$N), tolerance = 0))) {
    .archive_stop("The source theta or N disagrees with the selected configuration index.")
  }
  error <- pred_theta - theta
  se <- sqrt(pred_sig2 / N)
  critical <- stats::qnorm(.975)
  bias <- apply(error, 2L, stats::median)
  RMSE <- sqrt(colMeans(error^2))
  AL <- apply(2 * critical * se, 2L, stats::median)
  AC <- colMeans(abs(error) <= critical * se)
  ESD <- apply(pred_theta, 2L, function(x) 1.4826 * stats::median(abs(x - stats::median(x))))
  ASD <- apply(se, 2L, stats::median)
  sig2 <- apply(pred_sig2, 2L, stats::median)
  paper_method_names <- c("Oracle", "TR", "O-TR", "QR1", "QR2", "O-TR'", "MQR1", "MQR2")
  table_summary <- data.frame(method = method_names, bias = unname(bias),
    RMSE = unname(RMSE), AL = unname(AL), AC = unname(AC),
    ESD = unname(ESD), ASD = unname(ASD), row.names = NULL)
  table_summary <- table_summary[match(paper_method_names, table_summary$method), , drop = FALSE]
  rownames(table_summary) <- NULL
  results_par <- cbind(pred_theta, pred_sig2)
  colnames(results_par) <- c(paste0("pred_theta:", method_names), paste0("pred_sig2:", method_names))
  result <- list(source_id = as.character(case$id), source_path = as.character(case$path),
    group = case$group, config = case$config, replication_ids = replication_ids,
    tmax = 100L, N = N, theta = theta, seed = get("seed", envir = execution),
    num_totest = length(method_names), method_names = method_names,
    paper_method_names = paper_method_names, pred_theta = pred_theta,
    pred_sig2 = pred_sig2, results_par = results_par, sig2 = sig2,
    bias = bias, RMSE = RMSE, AL = AL, AC = AC, ESD = ESD, ASD = ASD,
    table_summary = table_summary)
  .archive_require_finite(result, "simulation result")
  result
}

.archive_real_result <- function(execution, replication_ids) {
  if (exists("method_names", envir = execution, inherits = FALSE)) {
    method_names <- get("method_names", envir = execution, inherits = FALSE)
  } else if (exists("pred_theta", envir = execution, inherits = FALSE)) {
    method_names <- colnames(get("pred_theta", envir = execution, inherits = FALSE))
  } else method_names <- NULL
  allowed <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR")
  if (!is.character(method_names) || !length(method_names) || anyNA(method_names) ||
      anyDuplicated(method_names) || any(!method_names %in% allowed)) {
    .archive_stop("Real-data output must name its retained methods explicitly; deprecated methods and Oracle cannot be exported.")
  }
  matrices <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
                "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
  result <- list(source_id = "real", source_path = "configurations/real_2.R",
    replication_ids = replication_ids, method_names = method_names,
    num_totest = length(method_names), tmax = 100L)
  for (name in matrices) {
    value <- .archive_matrix(execution, name, length(replication_ids), length(method_names))
    dimnames(value) <- list(as.character(replication_ids), method_names)
    if (startsWith(name, "pred_sig2") && any(value < 0, na.rm = TRUE)) {
      .archive_stop(name, " contains a negative influence-function variance.")
    }
    result[[name]] <- value
  }
  # Explicitly whitelist derived results. No patient-level data, functions, or
  # historical estimator objects enter the clean result.
  fields <- c("N_all", "N", "seed", "method_columns", "result_summary", "proportion",
    "sig2", "sig2_IDE", "sig2_DE", "sig2_TE", "theta", "theta0", "theta1",
    "DE", "IDE", "TE", "p_theta", "p_ide", "p_de", "p_te")
  for (name in fields) {
    if (exists(name, envir = execution, inherits = FALSE)) {
      value <- get(name, envir = execution, inherits = FALSE)
      if (is.function(value)) .archive_stop("Expected a result, but ", name, " is a function.")
      if (is.numeric(value) && is.null(dim(value)) && length(value) == length(method_names)) {
        names(value) <- method_names
      }
      result[[name]] <- value
    }
  }
  if (!is.data.frame(result$result_summary) || !"method" %in% names(result$result_summary) ||
      !setequal(as.character(result$result_summary$method), method_names) ||
      anyDuplicated(result$result_summary$method)) {
    .archive_stop("The real-data source must produce a named result_summary for every retained method.")
  }
  result$results_par <- do.call(cbind, result[matrices])
  colnames(result$results_par) <- unlist(lapply(matrices, function(metric) paste(metric, method_names, sep = ":")))
  result$paper_method_names <- c("TR", "O-TR", "QR1", "QR2", "O-TR'", "MQR1", "MQR2")
  result$paper_method_names <- result$paper_method_names[result$paper_method_names %in% method_names]
  result$table_summary <- result$result_summary[match(result$paper_method_names,
    as.character(result$result_summary$method)), , drop = FALSE]
  rownames(result$table_summary) <- NULL
  .archive_require_finite(result, "real-data result")
  result
}

run_archive <- function(group, config = NULL, replication_ids = 1:100,
                        archive_root = getwd()) {
  replication_ids <- .archive_replication_ids(replication_ids)
  archive_root <- normalizePath(archive_root, mustWork = TRUE)
  real <- identical(as.character(group), "real")
  if (real) {
    case <- list(id = "real", path = "configurations/real_2.R", tmax = 100L)
    output_name <- "real"
  } else {
    if (length(group) != 1L || !as.character(group) %in% as.character(1:6)) {
      .archive_stop("group must be one of 1 through 6, or 'real'.")
    }
    case <- .archive_select_case(archive_root, group, config)
    selector <- if (as.integer(group) <= 4L) "N" else if (as.integer(group) == 5L) "signal" else "s_level"
    output_name <- paste0("simu_", group, "_", selector, config)
  }
  source_path <- .archive_source_path(archive_root, case$path)
  output_root <- file.path(archive_root, "output", output_name)
  if (!dir.exists(output_root) && !dir.create(output_root, recursive = TRUE)) {
    .archive_stop("Cannot create output directory: ", output_root)
  }
  run_directory <- tempfile(pattern = paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_pid", Sys.getpid(), "_"), tmpdir = output_root)
  if (!dir.create(run_directory)) .archive_stop("Cannot create a fresh run directory.")
  run_id <- basename(run_directory)
  raw_directory <- file.path(run_directory, "internal_raw")
  if (!dir.create(raw_directory)) .archive_stop("Cannot create the raw-output directory.")
  if (real) {
    helper_path <- .archive_source_path(archive_root, "configurations/real_2_data_generate.R")
    if (!file.copy(helper_path, file.path(raw_directory, "real_2_data_generate.R"), overwrite = FALSE)) {
      .archive_stop("Cannot stage the real-data implementation helper.")
    }
  }
  manifest <- list(source_id = case$id, source_path = case$path,
    replication_ids = replication_ids, tmax = 100L, run_id = run_id,
    started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), status = "running")
  saveRDS(manifest, file.path(run_directory, "run_manifest.rds"))
  execution <- new.env(parent = .GlobalEnv)
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  on.exit({
    if (exists("cl", envir = execution, inherits = FALSE)) {
      cluster <- get("cl", envir = execution, inherits = FALSE)
      if (inherits(cluster, "cluster")) try(parallel::stopCluster(cluster), silent = TRUE)
    }
  }, add = TRUE)
  setwd(raw_directory)
  message("Run output: ", run_directory)
  result <- tryCatch({
    .archive_evaluate(source_path, replication_ids, execution, real = real)
    if (real) .archive_real_result(execution, replication_ids) else
      .archive_simulation_result(execution, case, replication_ids)
  }, error = function(error) {
    manifest$status <- "failed"
    manifest$error <- conditionMessage(error)
    saveRDS(manifest, file.path(run_directory, "run_manifest.rds"))
    .archive_stop("Archive run failed; generated diagnostics remain in ", raw_directory,
                  ".\n", conditionMessage(error))
  })
  result$run_id <- run_id
  id <- gsub("[^A-Za-z0-9_.-]", "_", as.character(case$id))
  final_path <- file.path(run_directory, paste0("Final_", id, "_", run_id, ".RData"))
  save(list = names(result), envir = list2env(result, parent = emptyenv()), file = final_path)
  summary_path <- file.path(run_directory, paste0("summary_", id, "_", run_id, ".csv"))
  utils::write.csv(result$table_summary, file = summary_path, row.names = FALSE)
  manifest$status <- "complete"
  manifest$finished_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  manifest$final_file <- basename(final_path)
  manifest$summary_file <- basename(summary_path)
  saveRDS(manifest, file.path(run_directory, "run_manifest.rds"))
  message("Clean result: ", final_path)
  invisible(list(final_file = final_path, summary_file = summary_path,
    run_directory = run_directory, manifest = manifest))
}

# Lightweight deterministic checks: no model fits, packages, or source data.
archive_runtime_self_test <- function() {
  execution <- new.env(parent = baseenv())
  execution$pred_theta <- matrix(seq_len(30L) / 10, nrow = 3L)
  execution$pred_sig2 <- matrix(4, nrow = 3L, ncol = 10L)
  execution$pred_theta[, c(4L, 9L)] <- NA_real_
  execution$pred_sig2[, c(4L, 9L)] <- NA_real_
  execution$theta <- 1
  execution$N <- 100
  execution$seed <- 1001
  case <- list(id = "synthetic", path = "synthetic.R", group = 1L,
               config = 100, theta = 1, N = 100, tmax = 100)
  result <- .archive_simulation_result(execution, case, c(2L, 7L, 100L))
  stopifnot(identical(dim(result$pred_theta), c(3L, 8L)),
    identical(result$replication_ids, c(2L, 7L, 100L)),
    identical(result$table_summary$method, result$paper_method_names),
    isTRUE(all.equal(unname(result$bias[1L]), -.8)),
    isTRUE(all.equal(unname(result$RMSE[1L]), sqrt(mean(c(-.9, -.8, -.7)^2)))),
    isTRUE(all.equal(unname(result$ASD[1L]), .2)),
    isTRUE(all.equal(unname(result$AL[1L]), 2 * stats::qnorm(.975) * .2)),
    isTRUE(all.equal(unname(result$ESD[1L]), 1.4826 * .1)))
  execution$pred_theta[1L, 1L] <- NA_real_
  rejected <- try(.archive_simulation_result(execution, case, c(2L, 7L, 100L)), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"))
  execution$replication_ids <- NULL
  eval(quote(rm(list = ls())), envir = execution)
  stopifnot(length(ls(execution)) == 0L)
  source <- textConnection(c("rm(list=ls())", "tmax <- 100", "replication_ids <- 1:100",
    "v <- tmax", "results_par <- seq_len(20)", "pred_theta <- results_par[,1:10]",
    "pred_sig2 <- results_par[,11:20]", "summary <- apply(results_par[,1:10],2,stats::median)"))
  on.exit(close(source), add = TRUE)
  .archive_evaluate(source, 7L, execution, real = TRUE)
  stopifnot(identical(execution$replication_ids, 7L), execution$tmax == 100,
    execution$v == 1, identical(dim(execution$pred_theta), c(1L, 10L)),
    identical(dim(execution$pred_sig2), c(1L, 10L)))
  stopifnot(!.archive_lhs_is(quote(parallel::clusterEvalQ(cl, 1)), "replication_ids"))
  real_execution <- new.env(parent = baseenv())
  real_execution$method_names <- c("MQR2", "TR")
  for (name in c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
                 "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")) {
    assign(name, matrix(seq_len(4L), nrow = 2L), envir = real_execution)
  }
  real_execution$theta <- c(1, 2)
  real_execution$result_summary <- data.frame(method = c("MQR2", "TR"), theta = c(1, 2))
  real_result <- .archive_real_result(real_execution, c(2L, 9L))
  stopifnot(identical(dim(real_result$results_par), c(2L, 14L)),
    identical(real_result$table_summary$method, c("TR", "MQR2")),
    identical(names(real_result$theta), c("MQR2", "TR")),
    !"Oracle" %in% real_result$method_names)
  invisible(TRUE)
}
