# Regenerate the seven numerical tables from the bundled replication-level results.
# This script uses base R only and does not fit nuisance models.
.table_root <- function() {
  paths <- vapply(sys.frames(), function(x) {
    if (is.character(x$ofile) && length(x$ofile) == 1L) x$ofile else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) {
    paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  }
  if (!length(paths)) stop("Run with source() or Rscript.")
  dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
}
archive_root <- .table_root()
source(file.path(archive_root, "tools", "verify_results.R"), local = TRUE)
verify_archive_results(archive_root)
table_directory <- file.path(archive_root, "tables")
dir.create(table_directory, showWarnings = FALSE)
load_result <- function(id) {
  env <- new.env(parent = baseenv())
  load(file.path(archive_root, "results", paste0(id, ".RData")), envir = env)
  as.list(env)
}
cases <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
table_stems <- c("main_table_5", "main_table_6", "main_table_7",
                 "supplement_table_4", "supplement_table_5", "supplement_table_6")
for (group in 1:6) {
  selected <- Filter(function(x) identical(as.integer(x$group), group), cases)
  summary_rows <- list()
  display_rows <- list()
  latex_rows <- character()
  for (case in selected) {
    result <- load_result(case$id)
    summary <- result$table_summary
    summary_rows[[case$id]] <- cbind(
      configuration = case$id, N = case$N, selector = case$config, summary)
    display <- data.frame(
      configuration = case$id, method = summary$method,
      Bias = sprintf("%.3f", summary$bias), RMSE = sprintf("%.3f", summary$RMSE),
      AL = sprintf("%.3f", summary$AL), AC = sprintf("%.2f", summary$AC),
      stringsAsFactors = FALSE)
    display_rows[[case$id]] <- display
    latex_rows <- c(latex_rows, paste0("% ", case$id),
      vapply(seq_len(nrow(display)), function(i) {
        paste0(paste(unlist(display[i, c("method", "Bias", "RMSE", "AL", "AC")]),
                     collapse = " & "), " \\\\")
      }, character(1)), "")
  }
  stem <- file.path(table_directory, table_stems[group])
  utils::write.csv(do.call(rbind, summary_rows), paste0(stem, ".csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, display_rows), paste0(stem, "_display.csv"), row.names = FALSE)
  writeLines(latex_rows, paste0(stem, "_rows.tex"))
}
real <- load_result("real")
aggregate <- recompute_real_results(real)
summary <- aggregate$result_summary[
  match(real$paper_method_names, aggregate$result_summary$method), , drop = FALSE]
utils::write.csv(summary, file.path(table_directory, "main_table_8.csv"), row.names = FALSE)
format_p <- function(x) ifelse(x < .001, sprintf("%.1e", x), sprintf("%.3f", x))
display <- data.frame(
  method = summary$method,
  theta = sprintf("%.3f", summary$theta), theta_length = sprintf("%.3f", summary$theta_length),
  NDE = sprintf("%.3f", summary$NDE), NDE_length = sprintf("%.3f", summary$NDE_length),
  NDE_p = format_p(summary$NDE_p), NIE = sprintf("%.3f", summary$NIE),
  NIE_length = sprintf("%.3f", summary$NIE_length), NIE_p = format_p(summary$NIE_p),
  proportion = sprintf("%.2f", 100 * summary$proportion), stringsAsFactors = FALSE)
# Check every printed string against the saved manuscript display, not only rounded numerics.
for (i in seq_len(nrow(real$manuscript_display))) {
  row <- real$manuscript_display[i, ]
  if (row$method == "Common ATE") next
  actual <- display[match(row$method, display$method), row$metric]
  stopifnot(identical(unname(actual), row$table_display))
}
utils::write.csv(display, file.path(table_directory, "main_table_8_display.csv"), row.names = FALSE)
latex_rows <- vapply(seq_len(nrow(display)), function(i) {
  row <- unlist(display[i, ])
  row[length(row)] <- paste0(row[length(row)], "\\%")
  paste0(paste(row, collapse = " & "), " \\\\")
}, character(1))
common <- data.frame(ATE = sprintf("%.3f", summary$ATE[1]),
  ATE_length = sprintf("%.3f", summary$ATE_length[1]), ATE_p = format_p(summary$ATE_p[1]))
utils::write.csv(common, file.path(table_directory, "main_table_8_common_ATE.csv"), row.names = FALSE)
writeLines(c(latex_rows, "", paste0("% Common ATE: ", paste(unlist(common), collapse = ", "))),
           file.path(table_directory, "main_table_8_rows.tex"))
message("Regenerated main Tables 5-8 and supplementary Tables 4-6 in: ", table_directory)
