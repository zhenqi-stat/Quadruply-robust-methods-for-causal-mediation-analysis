# Master workflow. Default: regenerate tables. --fit: rerun all 25 analyses.
.master_paths <- vapply(sys.frames(), function(x) {
  if (is.character(x$ofile) && length(x$ofile) == 1L) x$ofile else ""
}, character(1))
.master_paths <- .master_paths[nzchar(.master_paths)]
if (!length(.master_paths)) {
  .master_paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
}
if (!length(.master_paths)) stop("Run with source() or Rscript.")
archive_root <- dirname(normalizePath(tail(.master_paths, 1L), mustWork = TRUE))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(arguments, "--fit"))) stop("Usage: Rscript run_all.R [--fit]")
if (!"--fit" %in% arguments) {
  source(file.path(archive_root, "reproduce_tables.R"), local = TRUE)
} else {
  source(file.path(archive_root, "tools", "archive_runtime.R"), local = TRUE)
  cases <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
  for (case in cases) {
    message("Starting ", case$id, "; detailed progress is written to its log file.")
    run_archive(case$group, case$config, replication_ids = 1:100, archive_root = archive_root)
    gc(verbose = FALSE)
  }
  message("Starting ACTG175; detailed progress is written to its log file.")
  run_archive("real", replication_ids = 1:100, archive_root = archive_root)
}
