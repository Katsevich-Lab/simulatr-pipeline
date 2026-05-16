#!/usr/bin/env Rscript
# run_benchmark.R <study.rds> <method> <row_idx> <B_check> <B_in> \
#                 <max_gb> <max_hours> <benchmark_memory> \
#                 <R_SESSION_GB> <GB_PROPORTION_USED>
#
# Phase: benchmarking. Calls simulatr_benchmark() to get decomposed
# per-replicate memory and time numbers, then computes how many
# processors are needed to fit B replicates inside max_gb / max_hours.
# Emits two files:
#   benchmarking_info_<method>_<row_idx>.rds  — single-row tibble
#   proc_id_info_<method>_<row_idx>.csv       — one row per processor

suppressPackageStartupMessages(library(simulatr))
`%||%` <- function(a, b) if (is.null(a)) b else a

# Locate _simulatr_pipeline_helpers.R alongside this script (bin/ is on
# PATH inside a Nextflow task; the file lives next to the runner).
locate_helpers <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  fp <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(fp) > 0L && nzchar(fp[[1L]])) {
    cand <- file.path(dirname(normalizePath(fp[[1L]])),
                      "_simulatr_pipeline_helpers.R")
    if (file.exists(cand)) return(cand)
  }
  cand <- Sys.which("_simulatr_pipeline_helpers.R")[[1L]]
  if (nzchar(cand)) return(cand)
  stop("Cannot find _simulatr_pipeline_helpers.R on PATH or next to this script.")
}
source(locate_helpers())

args <- commandArgs(trailingOnly = TRUE)
study_fp           <- args[[1L]]
method             <- args[[2L]]
row_idx            <- as.integer(args[[3L]])
B_check            <- as.integer(args[[4L]])
B_in               <- as.integer(args[[5L]])
max_gb             <- as.numeric(args[[6L]])
max_hours          <- as.numeric(args[[7L]])
benchmark_memory   <- as.numeric(args[[8L]])
R_SESSION_GB       <- as.numeric(args[[9L]])
GB_PROPORTION_USED <- as.numeric(args[[10L]])

BYTES_PER_GB         <- 2 ^ 30
SECONDS_PER_HOUR     <- 60 * 60
HOURS_PROPORTION_USED <- 0.9

if (benchmark_memory == 0) max_gb <- Inf

with_error_sidecar("benchmark", row_idx, expr = {
  study <- readRDS(study_fp)
  one   <- study_with_one_method(study, method)

  bench <- simulatr_benchmark(one, row_idx,
                              n_warmup = 1L, n_measure = max(1L, B_check))
  pick <- function(phase) {
    row <- bench[bench$phase == phase, , drop = FALSE]
    if (nrow(row) == 0L) return(c(bytes = 0, seconds = 0))
    c(bytes = max(row$bytes), seconds = max(row$seconds))
  }
  gen <- pick("generate")
  run <- pick("run")
  ev  <- pick("evaluate")

  data_gb_per_rep   <- gen[["bytes"]]   / BYTES_PER_GB
  method_gb_per_rep <- run[["bytes"]]   / BYTES_PER_GB
  result_gb_per_rep <- ev[["bytes"]]    / BYTES_PER_GB
  if (result_gb_per_rep == 0) result_gb_per_rep <- max(method_gb_per_rep * 0.1, 1e-6)
  data_hours_per_rep   <- gen[["seconds"]] / SECONDS_PER_HOUR
  method_hours_per_rep <- run[["seconds"]] / SECONDS_PER_HOUR

  B <- if (B_in != 0) B_in else study@fixed_parameters$B

  max_reps_gb <- (max_gb * GB_PROPORTION_USED -
                  R_SESSION_GB - data_gb_per_rep - method_gb_per_rep) /
                 max(result_gb_per_rep, 1e-9)
  max_reps_hours <- max_hours * HOURS_PROPORTION_USED /
                    max(data_hours_per_rep + method_hours_per_rep, 1e-9)

  if (max_reps_gb < 1) {
    stop(sprintf(
      "One rep of method '%s' (grid %d) requires %.3f GB but max_gb is %.2f.",
      method, row_idx,
      R_SESSION_GB + data_gb_per_rep + method_gb_per_rep + result_gb_per_rep,
      max_gb))
  }
  if (max_reps_hours < 1) {
    stop(sprintf(
      "One rep of method '%s' (grid %d) requires %.3f hours but max_hours is %.2f.",
      method, row_idx, data_hours_per_rep + method_hours_per_rep, max_hours))
  }

  max_reps     <- min(max_reps_gb, max_reps_hours)
  n_processors <- max(1L, ceiling(B / max_reps))

  benchmarking_info <- data.frame(
    method               = method,
    grid_id              = row_idx,
    data_gb_per_rep      = data_gb_per_rep,
    method_gb_per_rep    = method_gb_per_rep,
    result_gb_per_rep    = result_gb_per_rep,
    method_hours_per_rep = method_hours_per_rep,
    data_hours_per_rep   = data_hours_per_rep,
    max_reps_gb          = max_reps_gb,
    max_reps_hours       = max_reps_hours,
    max_reps             = max_reps,
    n_processors         = n_processors
  )
  saveRDS(benchmarking_info,
          sprintf("benchmarking_info_%s_%d.rds", method, row_idx))

  proc_id_info <- data.frame(method = method,
                             grid_id = row_idx,
                             proc_id = seq_len(n_processors),
                             n_processors = n_processors)
  write.table(proc_id_info,
              file = sprintf("proc_id_info_%s_%d.csv", method, row_idx),
              col.names = FALSE, row.names = FALSE,
              quote = FALSE, sep = ",")
})
