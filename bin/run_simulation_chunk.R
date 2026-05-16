#!/usr/bin/env Rscript
# run_simulation_chunk.R <study.rds> <method> <row_idx> \
#                       <proc_id> <n_processors> <B_in>
#
# Phase: simulation. Runs the slice of replicates assigned to this
# processor for a single (method, grid_row). Mirrors simulatr_run()'s
# interleaved generate-then-run-then-discard loop, restricted to the
# proc_id_b subset.
#
# Emits chunk_result_<method>_<row_idx>_<proc_id>.rds with the
# {grid_id, method, chunk_id, run_id, output} tibble shape from the
# package's `results` schema.

suppressPackageStartupMessages(library(simulatr))
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
proc_id      <- as.integer(args[[4L]])
n_processors <- as.integer(args[[5L]])
B_in         <- as.integer(args[[6L]])

with_error_sidecar("chunk", row_idx, chunk_id = proc_id, expr = {
  study <- readRDS(study_fp)

  if (!(method_name %in% names(study@methods))) {
    stop(sprintf("Method '%s' not in study (have: %s).",
                 method_name, paste(names(study@methods), collapse = ", ")))
  }
  method <- study@methods[[method_name]]
  model  <- study@models[[1L]]
  if (model@batched || method@batched) {
    stop("batched = TRUE not supported by run_simulation_chunk.R yet ",
         "(use the single-process simulatr_run on the laptop, or set ",
         "n_processors = 1 and rerun).")
  }

  B    <- if (B_in != 0L) B_in else study@fixed_parameters$B
  seed <- study@fixed_parameters$seed
  all_b     <- seq_len(B)
  proc_id_b <- all_b[((all_b - 1L) %% n_processors) + 1L == proc_id]

  ground_truth <- if ("ground_truth" %in% colnames(study@parameter_grid)) {
    study@parameter_grid[["ground_truth"]][[row_idx]]
  } else NULL

  pull_arg <- function(name) {
    if (name %in% colnames(study@parameter_grid)) {
      study@parameter_grid[[name]][[row_idx]]
    } else if (name %in% names(study@fixed_parameters)) {
      study@fixed_parameters[[name]]
    } else if (name == "ground_truth") {
      ground_truth
    } else {
      stop(sprintf("Cannot resolve formal '%s' against parameter_grid or fixed_parameters.", name))
    }
  }
  build_args <- function(fn, skip_first = FALSE) {
    fl <- names(formals(fn))
    if (skip_first) fl <- fl[-1L]
    out <- vector("list", length(fl))
    names(out) <- fl
    for (a in fl) out[[a]] <- pull_arg(a)
    out
  }

  for (p in unique(c(model@packages, method@packages))) {
    if (nzchar(p)) requireNamespace(p, quietly = TRUE)
  }

  args_gen <- build_args(model@generate)
  args_run <- build_args(method@run, skip_first = TRUE)

  result_list <- vector("list", length(proc_id_b))
  for (k in seq_along(proc_id_b)) {
    b <- proc_id_b[[k]]
    draw <- R.utils::withSeed(do.call(model@generate, args_gen),
                              seed = seed + b)
    captured <- tryCatch(
      withCallingHandlers(
        R.utils::withSeed(do.call(method@run, c(list(draw), args_run)),
                          seed = seed + b),
        warning = function(w) {
          message(sprintf("(b=%d) warning captured: %s", b, conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        write_error_json("chunk", row_idx, chunk_id = proc_id,
                         run_id = b, e)
        structure(list(message = conditionMessage(e),
                       call = paste(deparse(conditionCall(e)),
                                    collapse = " ")),
                  class = c("simulatr_error", "list"))
      }
    )
    result_list[[k]] <- tibble::tibble(
      grid_id  = as.integer(row_idx),
      method   = method_name,
      chunk_id = as.integer(proc_id),
      run_id   = as.integer(b),
      output   = list(captured)
    )
    rm(draw); gc(verbose = FALSE)
  }

  out_df <- dplyr::bind_rows(result_list)
  out_fp <- sprintf("chunk_result_%s_%d_%d.rds", method_name, row_idx, proc_id)
  saveRDS(out_df, out_fp)
})
