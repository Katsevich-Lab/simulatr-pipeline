#!/usr/bin/env Rscript
# collect_memory.R <out_tsv> [<out_parquet>]
#
# Concatenates every .command.memory.tsv staged in CWD by the
# Layer-2 hooks (Slurm `sacct` or rss_sampler.sh) into one TSV /
# parquet, normalising the rows into the contract's `memory` schema
# (grid_id, method, kind, phase, metric, value).

args <- commandArgs(trailingOnly = TRUE)
out_tsv     <- args[[1L]]
out_parquet <- if (length(args) >= 2L) args[[2L]] else NA_character_

files <- list.files(pattern = "\\.command\\.memory\\.tsv$", recursive = TRUE,
                    all.files = TRUE, full.names = TRUE)
if (length(files) == 0L) {
  cat("No .command.memory.tsv files found; emitting empty TSV.\n")
  writeLines("grid_id\tmethod\tkind\tphase\tmetric\tvalue", out_tsv)
  quit(status = 0)
}

parse_kb <- function(s) {
  s <- gsub("[KMG]$", "", s, perl = TRUE)
  suppressWarnings(as.numeric(s))
}

rows <- lapply(files, function(f) {
  txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  if (length(txt) < 2L) return(NULL)
  hdr <- strsplit(txt[[1L]], "[|\t]")[[1L]]
  vals <- strsplit(txt[[2L]], "[|\t]")[[1L]]
  d <- setNames(as.list(vals), hdr)
  raw <- d$MaxRSS
  if (is.null(raw) || !nzchar(raw)) return(NULL)
  unit_mult <- if (grepl("K$", raw)) 1024 else if (grepl("M$", raw)) 1024^2 else if (grepl("G$", raw)) 1024^3 else 1
  bytes <- parse_kb(raw) * unit_mult
  if (is.na(bytes)) return(NULL)
  data.frame(
    grid_id = NA_integer_,
    method  = NA_character_,
    kind    = if (grepl("local$", d$JobName %||% "")) "rss_sampler" else "sacct",
    phase   = NA_character_,
    metric  = "max_rss_gb",
    value   = bytes / 2 ^ 30,
    stringsAsFactors = FALSE
  )
})
`%||%` <- function(a, b) if (is.null(a)) b else a
rows <- Filter(Negate(is.null), rows)

if (length(rows) == 0L) {
  writeLines("grid_id\tmethod\tkind\tphase\tmetric\tvalue", out_tsv)
  quit(status = 0)
}

mem <- do.call(rbind, rows)
write.table(mem, file = out_tsv, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

if (!is.na(out_parquet) && requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(mem, out_parquet)
}
