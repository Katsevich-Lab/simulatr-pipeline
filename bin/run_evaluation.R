#!/usr/bin/env Rscript
# run_evaluation.R <study.rds> <method> <row_idx> \
#                  <chunk_results_glob...> <benchmarking_info.rds>
#
# Phase: evaluation. Concatenates this (method, row_idx)'s chunk
# results, joins the benchmarking info, applies the study's metrics
# via simulatr_evaluate(), and writes
# <method>_<row_idx>_results.rds = list(results, metrics).

suppressPackageStartupMessages({
  library(simulatr)
  library(dplyr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

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
study_fp     <- args[[1L]]
method_name  <- args[[2L]]
row_idx      <- as.integer(args[[3L]])
nextflow_fps <- args[4L:length(args)]

chunk_fps <- nextflow_fps[grepl("chunk_result", nextflow_fps)]
bench_fp  <- nextflow_fps[grepl("benchmarking_info", nextflow_fps)][[1L]]

with_error_sidecar("evaluate", row_idx, expr = {
  study <- readRDS(study_fp)

  results <- lapply(chunk_fps, readRDS) |>
    data.table::rbindlist() |>
    dplyr::as_tibble()

  benchmarking_info <- readRDS(bench_fp) |> dplyr::as_tibble()

  evals <- simulatr_evaluate(study, results)
  summary_tbl <- simulatr_evals_summarise(evals)

  bench_long <- benchmarking_info |>
    dplyr::rename(gb_per_rep  = method_gb_per_rep,
                  hrs_per_rep = method_hours_per_rep) |>
    dplyr::select(method, grid_id, gb_per_rep, hrs_per_rep, n_processors) |>
    tidyr::pivot_longer(c(gb_per_rep, hrs_per_rep, n_processors),
                        names_to  = "metric",
                        values_to = "mean")

  grid_no_gt <- study@parameter_grid |>
    dplyr::mutate(grid_id = dplyr::row_number())
  if ("ground_truth" %in% colnames(grid_no_gt)) {
    grid_no_gt <- dplyr::select(grid_no_gt, -ground_truth)
  }

  if (nrow(summary_tbl) > 0L) {
    metrics <- summary_tbl |>
      dplyr::transmute(method, grid_id, metric, mean, se) |>
      dplyr::bind_rows(bench_long) |>
      dplyr::left_join(grid_no_gt, by = "grid_id") |>
      dplyr::select(-grid_id) |>
      dplyr::relocate(method, metric, mean, se)
  } else {
    metrics <- bench_long |>
      dplyr::left_join(grid_no_gt, by = "grid_id") |>
      dplyr::select(-grid_id) |>
      dplyr::relocate(method, metric, mean)
  }

  output <- list(
    results = results |> dplyr::select(method, grid_id, run_id, output),
    metrics = metrics
  )
  saveRDS(output, sprintf("%s_%d_results.rds", method_name, row_idx))
})
