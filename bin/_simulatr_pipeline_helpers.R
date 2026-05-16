# Shared helpers used by the bin/*.R pipeline scripts.
#
# Sourced via `source("/path/to/_simulatr_pipeline_helpers.R")` from
# each bin/*.R script (Nextflow puts bin/ on PATH, so the path is
# resolvable as Sys.which("_simulatr_pipeline_helpers.R")).

# Convert any captured condition into a structured JSON sidecar.
write_error_json <- function(phase, grid_id, chunk_id = NA_integer_,
                             run_id = NA_integer_, cond) {
  payload <- list(
    phase     = phase,
    grid_id   = grid_id,
    chunk_id  = chunk_id,
    run_id    = run_id,
    type      = paste(class(cond), collapse = "/"),
    message   = conditionMessage(cond),
    traceback = paste(utils::limitedLabels(sys.calls()), collapse = "\n")
  )
  out <- file.path(
    getwd(),
    sprintf("error_%s_%s_%s.json", phase, grid_id,
            if (is.na(chunk_id)) "0" else chunk_id)
  )
  writeLines(
    jsonify(payload),
    out
  )
  invisible(out)
}

# Tiny dependency-free JSON serializer for the small payloads we emit.
jsonify <- function(x) {
  esc <- function(s) {
    s <- gsub("\\\\", "\\\\\\\\", as.character(s))
    s <- gsub("\"",   "\\\\\"",   s)
    s <- gsub("\n",   "\\\\n",    s)
    paste0("\"", s, "\"")
  }
  format_one <- function(k, v) {
    if (is.null(v) || (length(v) == 1 && is.na(v))) {
      paste0(esc(k), ": null")
    } else if (is.numeric(v) && length(v) == 1) {
      paste0(esc(k), ": ", v)
    } else {
      paste0(esc(k), ": ", esc(v))
    }
  }
  body <- vapply(names(x), function(k) format_one(k, x[[k]]), character(1))
  paste0("{", paste(body, collapse = ", "), "}")
}

# Run an expression, write error.json on any error, then re-throw so
# Nextflow sees a non-zero exit.
with_error_sidecar <- function(phase, grid_id, chunk_id = NA, run_id = NA, expr) {
  expr <- substitute(expr)
  env  <- parent.frame()
  tryCatch(
    eval(expr, env),
    error = function(e) {
      write_error_json(phase, grid_id, chunk_id, run_id, e)
      stop(e)
    }
  )
}

# Build a simulatr_study with a single named method for sub-runs that
# only operate on one method at a time. Reused by run_benchmark.R,
# run_simulation_chunk.R, run_evaluation.R.
study_with_one_method <- function(study, method_name) {
  if (!(method_name %in% names(study@methods))) {
    stop(sprintf("Method '%s' not found in study (have: %s).",
                 method_name,
                 paste(names(study@methods), collapse = ", ")))
  }
  # simulatr_study's S4 validity requires plain data.frame (not tbl_df).
  grid <- as.data.frame(study@parameter_grid, check.names = FALSE)
  simulatr::simulatr_study(
    parameter_grid   = grid,
    fixed_parameters = study@fixed_parameters,
    models           = study@models,
    methods          = study@methods[method_name],
    metrics          = study@metrics
  )
}
