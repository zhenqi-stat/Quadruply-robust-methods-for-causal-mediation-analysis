# Check code and saved results by default; add --dgp to regenerate the 24 DGP pools.
.verification_root <- function() {
  frames <- sys.frames()
  paths <- vapply(frames, function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (length(paths)) {
    path <- paths[[length(paths)]]
  } else {
    path <- commandArgs(trailingOnly = FALSE)
    path <- path[startsWith(path, "--file=")]
    if (!length(path)) stop("Run verify.R with Rscript or source().")
    path <- sub("^--file=", "", path[[1L]])
  }
  if (!file.exists(path) && file.exists(basename(path))) path <- basename(path)
  dirname(normalizePath(path, mustWork = TRUE))
}

.has_constant_true_if <- function(expression) {
  if (!is.call(expression)) return(FALSE)
  if (identical(expression[[1L]], as.name("if")) &&
      identical(expression[[2L]], TRUE)) return(TRUE)
  any(vapply(as.list(expression), .has_constant_true_if, logical(1)))
}

archive_root <- .verification_root()
checksum_file <- file.path(archive_root, "file_checksums.csv")
checksums <- utils::read.csv(checksum_file, stringsAsFactors = FALSE)
stopifnot(nrow(checksums) > 25L, !anyDuplicated(checksums$path))
paths <- file.path(archive_root, checksums$path)
if (any(!file.exists(paths))) stop("An archived file is missing.")
actual <- unname(tools::md5sum(paths))
if (!identical(actual, checksums$md5)) {
  stop("Checksum mismatch: ", paste(checksums$path[actual != checksums$md5], collapse = ", "))
}

source(file.path(archive_root, "tools", "archive_runtime.R"), local = TRUE)
archive_runtime_self_test()
index <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
stopifnot(length(index) == 24L)
for (case in index) {
  selected <- .archive_select_case(archive_root, case$group, case$config)
  stopifnot(identical(selected$id, case$id), case$tmax == 100)
  path <- .archive_source_path(archive_root, case$path)
  expressions <- parse(path, keep.source = FALSE)
  stopifnot(!any(vapply(expressions, .has_constant_true_if, logical(1))))
  active_text <- paste(deparse(expressions, width.cutoff = 500L), collapse = "\n")
  completions <- c("MQR1 accomplished", "MQR2 accomplished",
    "QR1 accomplished", "QR2 accomplished", "O-TR' accomplished",
    "TR & O-TR accomplished")
  stopifnot(all(vapply(completions, grepl, logical(1), x = active_text, fixed = TRUE)))
  stopifnot(sum(vapply(expressions, .archive_lhs_is, logical(1), name = "replication_ids")) == 1L)
}
r_files <- list.files(archive_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("/output/", r_files, fixed = TRUE)]
for (path in r_files) parse(path, keep.source = FALSE)
message("PASS: checksums, ", length(r_files), " R source files, 24 configurations, method switches, and runtime self-tests.")

source(file.path(archive_root, "tools", "verify_results.R"), local = TRUE)
verify_archive_results(archive_root)

if ("--dgp" %in% commandArgs(trailingOnly = TRUE)) {
  source(file.path(archive_root, "tools", "verify_dgp.R"), local = TRUE)
  result <- verify_archive_dgp(archive_root)
  print(result, row.names = FALSE)
  message("PASS: 24 finite-pool targets and 2,400 Oracle estimates/variances; no nuisance-model fitting.")
}
