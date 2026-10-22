# Simulation 5. Change only the selector and requested replication IDs below.
signal <- 1
# This is an execution subset; the scientific DGP always retains tmax = 100.
replication_ids <- 1:100

.archive_entry_root <- function() {
  frames <- sys.frames()
  source_files <- vapply(frames, function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files)) {
    path <- source_files[[length(source_files)]]
  } else {
    file_argument <- commandArgs(trailingOnly = FALSE)
    file_argument <- file_argument[startsWith(file_argument, "--file=")]
    if (!length(file_argument)) stop("Run this file with Rscript or source().", call. = FALSE)
    path <- sub("^--file=", "", file_argument[[1L]])
  }
  if (!file.exists(path) && file.exists(basename(path))) path <- basename(path)
  dirname(normalizePath(path, mustWork = TRUE))
}
archive_root <- .archive_entry_root()
source(file.path(archive_root, "tools", "archive_runtime.R"), local = TRUE)
run_archive(group = 5L, config = signal, replication_ids = replication_ids,
            archive_root = archive_root)

