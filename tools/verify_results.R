# Recompute the ACTG175 table from the seven saved split-level matrices.
# No patient-level data or nuisance estimation is required.
recompute_real_results <- function(record) {
  metrics <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
               "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
  if (!is.list(record) || !all(c("N", metrics) %in% names(record))) {
    stop("The real-data result must contain N and all seven result matrices.")
  }
  if (length(record$N) != 1L || !is.finite(record$N) || record$N <= 0) {
    stop("N must be a positive finite scalar.")
  }
  reference_dim <- dim(record$pred_theta)
  methods <- colnames(record$pred_theta)
  if (length(reference_dim) != 2L || reference_dim[1] != 100L ||
      reference_dim[2] < 1L || is.null(methods) || anyDuplicated(methods)) {
    stop("The result matrices must contain 100 splits and named method columns.")
  }
  for (metric in metrics) {
    x <- record[[metric]]
    if (!is.matrix(x) || !is.numeric(x) ||
        !identical(dim(x), reference_dim) ||
        !identical(colnames(x), methods) || !all(is.finite(x))) {
      stop("Invalid or incomplete result matrix: ", metric)
    }
    if (grepl("^pred_sig2", metric) && any(x < 0)) {
      stop("Negative variance in ", metric)
    }
  }
  out <- list()
  out$DE_split <- record$pred_theta - record$pred_theta0
  out$IDE_split <- record$pred_theta1 - record$pred_theta
  out$TE_split <- record$pred_theta1 - record$pred_theta0
  for (mapping in list(c("theta", "pred_theta"), c("theta0", "pred_theta0"),
                       c("theta1", "pred_theta1"))) {
    out[[mapping[1]]] <- apply(record[[mapping[2]]], 2, median)
  }
  for (mapping in list(c("DE", "DE_split"), c("IDE", "IDE_split"), c("TE", "TE_split"))) {
    out[[mapping[1]]] <- apply(out[[mapping[2]]], 2, median)
  }
  # Apply the saved analysis's scalar split-aggregation formula unchanged.
  out$sig2 <- apply(record$pred_sig2 + sweep(record$pred_theta, 2, out$theta, "-")^2, 2, median)
  out$sig2_DE <- apply(record$pred_sig2_DE + sweep(out$DE_split, 2, out$DE, "-")^2, 2, median)
  out$sig2_IDE <- apply(record$pred_sig2_IDE + sweep(out$IDE_split, 2, out$IDE, "-")^2, 2, median)
  out$sig2_TE <- apply(record$pred_sig2_TE + sweep(out$TE_split, 2, out$TE, "-")^2, 2, median)
  out$p_theta <- 2 * pnorm(-abs(out$theta / sqrt(out$sig2 / record$N)))
  out$p_de <- 2 * pnorm(-abs(out$DE / sqrt(out$sig2_DE / record$N)))
  out$p_ide <- 2 * pnorm(-abs(out$IDE / sqrt(out$sig2_IDE / record$N)))
  out$p_te <- 2 * pnorm(-abs(out$TE / sqrt(out$sig2_TE / record$N)))
  out$proportion <- out$IDE / out$TE
  out$result_summary <- data.frame(
    method = methods,
    theta = unname(out$theta),
    theta_length = 2 * qnorm(.975) * sqrt(out$sig2 / record$N),
    NDE = unname(out$DE),
    NDE_length = 2 * qnorm(.975) * sqrt(out$sig2_DE / record$N),
    NDE_p = unname(out$p_de),
    NIE = unname(out$IDE),
    NIE_length = 2 * qnorm(.975) * sqrt(out$sig2_IDE / record$N),
    NIE_p = unname(out$p_ide),
    proportion = unname(out$proportion),
    ATE = unname(out$TE),
    ATE_length = 2 * qnorm(.975) * sqrt(out$sig2_TE / record$N),
    ATE_p = unname(out$p_te),
    row.names = methods
  )
  out$confidence_intervals <- lapply(list(
    theta10 = c("theta", "sig2"), NDE = c("DE", "sig2_DE"),
    NIE = c("IDE", "sig2_IDE"), ATE = c("TE", "sig2_TE")
  ), function(mapping) {
    half_length <- qnorm(.975) * sqrt(out[[mapping[2]]] / record$N)
    cbind(lower = out[[mapping[1]]] - half_length,
          upper = out[[mapping[1]]] + half_length)
  })
  out
}

# Check the bundled final tables without fitting any nuisance model.
verify_archive_results <- function(archive_root) {
  archive_root <- normalizePath(archive_root, mustWork = TRUE)
  directory <- file.path(archive_root, "results")
  cases <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
  expected_files <- c(paste0(names(cases), ".RData"), "real.RData")
  files <- list.files(directory, pattern = "[.]RData$")
  stopifnot(length(files) == 25L, setequal(files, expected_files))
  methods <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
  equal <- function(x, y) {
    stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12)))
  }
  loaded <- function(name) {
    env <- new.env(parent = baseenv())
    load(file.path(directory, name), envir = env)
    as.list(env)
  }
  checked <- 0L
  for (id in names(cases)) {
    result <- loaded(paste0(id, ".RData"))
    case <- cases[[id]]
    stopifnot(identical(result$source_id, id),
              identical(result$replication_ids, 1:100),
              identical(result$method_names, methods),
              identical(dim(result$pred_theta), c(100L, 8L)),
              identical(dim(result$pred_sig2), c(100L, 8L)),
              identical(colnames(result$pred_theta), methods),
              identical(colnames(result$pred_sig2), methods),
              identical(rownames(result$pred_theta), as.character(1:100)),
              identical(rownames(result$pred_sig2), as.character(1:100)),
              all(is.finite(result$pred_theta)), all(is.finite(result$pred_sig2)),
              all(result$pred_sig2 >= 0), result$N == case$N,
              result$seed == case$seed, result$tmax == 100L)
    equal(result$theta, case$theta)
    error <- result$pred_theta - result$theta
    se <- sqrt(result$pred_sig2 / result$N)
    summary <- data.frame(
      method = methods, bias = unname(apply(error, 2, median)),
      RMSE = unname(sqrt(colMeans(error^2))),
      AL = unname(apply(2 * qnorm(.975) * se, 2, median)),
      AC = unname(colMeans(abs(error) <= qnorm(.975) * se)),
      ESD = unname(apply(result$pred_theta, 2,
        function(x) 1.4826 * median(abs(x - median(x))))),
      ASD = unname(apply(se, 2, median)), row.names = NULL)
    summary <- summary[match(result$paper_method_names, summary$method), , drop = FALSE]
    rownames(summary) <- NULL
    equal(result$table_summary, summary)
    for (field in c("bias", "RMSE", "AL", "AC", "ESD", "ASD")) {
      equal(unname(result[[field]]), summary[[field]][match(methods, summary$method)])
    }
    equal(unname(result$sig2), unname(apply(result$pred_sig2, 2, median)))
    equal(unname(result$results_par), unname(cbind(result$pred_theta, result$pred_sig2)))
    stopifnot(identical(dim(result$results_par), c(100L, 16L)),
              identical(result$manuscript_table$method, summary$method))
    for (field in c("bias", "RMSE", "AL", "AC")) {
      digits <- if (field == "AC") 2L else 3L
      stopifnot(all(round(summary[[field]], digits) == result$manuscript_table[[field]]))
    }
    map <- result$source_map
    stopifnot(nrow(map) == 800L,
              !anyDuplicated(paste(map$replication_id, map$method)),
              setequal(map$method, methods), all(table(map$method) == 100L),
              all(map$replication_id %in% 1:100), all(map$source_row >= 1L),
              all(map$source_file == basename(map$source_file)),
              all(grepl("^[0-9a-f]{32}$", map$source_md5)))
    checked <- checked + 32L
  }
  real <- loaded("real.RData")
  stopifnot(identical(real$source_id, "real"),
            identical(real$replication_ids, 1:100),
            identical(real$method_names, c("MQR2", "O-TR'", "QR2", "TR", "O-TR")),
            identical(dim(real$results_par), c(100L, 35L)))
  recomputed <- recompute_real_results(real)
  for (field in names(recomputed)) equal(real[[field]], recomputed[[field]])
  matrices <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
                "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
  equal(unname(real$results_par), unname(do.call(cbind, real[matrices])))
  equal(real$table_summary, {
    x <- recomputed$result_summary[match(real$paper_method_names,
      recomputed$result_summary$method), , drop = FALSE]
    rownames(x) <- NULL
    x
  })
  recipe <- real$source_recipe
  stopifnot(nrow(recipe) == 5L, identical(recipe$method, real$method_names),
            all(recipe$source_row_start == 1L), all(recipe$source_row_end == 100L),
            all(recipe$replication_id_start == 1L), all(recipe$replication_id_end == 100L),
            all(recipe$source_file == basename(recipe$source_file)),
            all(grepl("^[0-9a-f]{32}$", recipe$source_md5)))
  display <- real$manuscript_display
  stopifnot(nrow(display) == 48L, is.character(display$table_display))
  for (i in seq_len(nrow(display))) {
    row <- display[i, ]
    method <- if (row$method == "Common ATE") real$method_names[[1L]] else row$method
    value <- recomputed$result_summary[method, row$metric]
    if (row$metric == "proportion") value <- 100 * value
    expected <- row$table_display
    scientific <- grepl("e", expected, fixed = TRUE)
    mantissa <- if (scientific) strsplit(expected, "e", fixed = TRUE)[[1L]][[1L]] else expected
    fractional <- strsplit(mantissa, ".", fixed = TRUE)[[1L]]
    digits <- if (length(fractional) > 1L) nchar(fractional[[2L]]) else 0L
    actual <- sprintf(paste0("%.", digits, if (scientific) "e" else "f"), value)
    stopifnot(identical(actual, expected))
  }
  checked <- checked + nrow(display)
  message("PASS: 25 final RData files, 100 repetitions/splits each, and ",
          checked, " displayed table values recomputed.")
  invisible(TRUE)
}
