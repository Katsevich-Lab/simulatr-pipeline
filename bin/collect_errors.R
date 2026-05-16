#!/usr/bin/env Rscript
# collect_errors.R <out_tsv> [<out_parquet>]
#
# Concatenates every error_*.json staged in CWD by the with_error_sidecar
# helper into one TSV / parquet matching the contract's `errors` schema
# (phase, grid_id, chunk_id, run_id, type, message, traceback).

args <- commandArgs(trailingOnly = TRUE)
out_tsv     <- args[[1L]]
out_parquet <- if (length(args) >= 2L) args[[2L]] else NA_character_

files <- list.files(pattern = "^error_.*\\.json$", recursive = TRUE,
                    all.files = TRUE, full.names = TRUE)
if (length(files) == 0L) {
  cat("No error_*.json sidecars found; emitting empty TSV.\n")
  writeLines("phase\tgrid_id\tchunk_id\trun_id\ttype\tmessage\ttraceback",
             out_tsv)
  quit(status = 0)
}

# Tiny dependency-free JSON parser that handles the flat shape we emit.
parse_simple_json <- function(s) {
  s <- trimws(s)
  s <- sub("^\\{", "", s); s <- sub("\\}$", "", s)
  pairs <- regmatches(s, gregexpr('"[^"]+"\\s*:\\s*("[^"]*"|null|-?[0-9.]+)', s))[[1L]]
  out <- list()
  for (p in pairs) {
    parts <- regmatches(p, regexec('^"([^"]+)"\\s*:\\s*(.*)$', p))[[1L]]
    if (length(parts) < 3L) next
    key <- parts[[2L]]
    val_raw <- parts[[3L]]
    val <- if (val_raw == "null") NA
           else if (substr(val_raw, 1, 1) == "\"")
             gsub("\\\\n", "\n", gsub("\\\\\"", "\"", substr(val_raw, 2, nchar(val_raw) - 1)))
           else suppressWarnings(as.numeric(val_raw))
    out[[key]] <- val
  }
  out
}

rows <- lapply(files, function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  rec <- parse_simple_json(txt)
  data.frame(
    phase     = rec$phase     %||% NA_character_,
    grid_id   = as.integer(rec$grid_id  %||% NA),
    chunk_id  = as.integer(rec$chunk_id %||% NA),
    run_id    = as.integer(rec$run_id   %||% NA),
    type      = rec$type      %||% NA_character_,
    message   = rec$message   %||% NA_character_,
    traceback = rec$traceback %||% NA_character_,
    stringsAsFactors = FALSE
  )
})
`%||%` <- function(a, b) if (is.null(a)) b else a

errs <- do.call(rbind, rows)
write.table(errs, file = out_tsv, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

if (!is.na(out_parquet) && requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(errs, out_parquet)
}
