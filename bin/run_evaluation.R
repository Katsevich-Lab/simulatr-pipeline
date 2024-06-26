#!/usr/bin/env Rscript
library(simulatr)

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
simulatr_spec <- readRDS(args[1])
result_file_name <- args[2]
nextflow_res_filenames <- args[3:length(args)]

# extract command line arguments corresponding to simulation chunk results
chunk_result_filenames <- nextflow_res_filenames[grepl("chunk_result", nextflow_res_filenames)]

# extract command line arguments corresponding to benchmarking info
benchmarking_info_filenames <- nextflow_res_filenames[grepl("benchmarking_info", nextflow_res_filenames)]

# read the simulation chunk results
results <- lapply(X = chunk_result_filenames, FUN = readRDS) |>
  data.table::rbindlist() |>
  dplyr::as_tibble()

# read the benchmarking info
benchmarking_info <- lapply(X = benchmarking_info_filenames, FUN = readRDS) |>
  data.table::rbindlist() |>
  dplyr::as_tibble()

# join the results with the parameter grid
results_joined <- results |>
  dplyr::left_join(simulatr_spec@parameter_grid |>
    dplyr::mutate(grid_id = dplyr::row_number()) |>
    dplyr::select(grid_id, ground_truth),
  by = "grid_id"
  )

# evaluate the metrics
if (length(simulatr_spec@evaluation_functions) > 0) {
  metrics <- lapply(names(simulatr_spec@evaluation_functions), function(fun_name) {
    results_joined |>
      dplyr::rowwise() |>
      dplyr::mutate(metric = fun_name, value = simulatr_spec@evaluation_functions[[fun_name]](output, ground_truth)) |>
      dplyr::ungroup()
  }) |>
    data.table::rbindlist() |>
    dplyr::group_by(grid_id, method, metric) |>
    dplyr::summarise(mean = mean(value), se = sd(value) / sqrt(dplyr::n()), .groups = "drop") |>
    dplyr::bind_rows(benchmarking_info |>
      tidyr::pivot_longer(c(gb_per_rep, hrs_per_rep, n_processors),
        names_to = "metric",
        values_to = "mean"
      )) |>
    dplyr::left_join(simulatr_spec@parameter_grid |>
      dplyr::mutate(grid_id = dplyr::row_number()) |>
      dplyr::select(-ground_truth),
    by = "grid_id"
    ) |>
    dplyr::select(-grid_id) |>
    dplyr::relocate(method, metric, mean, se)
} else {
  metrics <- benchmarking_info |>
    tidyr::pivot_longer(c(gb_per_rep, hrs_per_rep, n_processors),
      names_to = "metric", values_to = "mean"
    ) |>
    dplyr::left_join(simulatr_spec@parameter_grid |>
      dplyr::mutate(grid_id = dplyr::row_number()) |>
      dplyr::select(-ground_truth),
    by = "grid_id"
    ) |>
    dplyr::select(-grid_id) |>
    dplyr::relocate(method, metric, mean)
}

# return
output <- list(
  results = results |> dplyr::select(method, grid_id, run_id, output),
  metrics = metrics
)
saveRDS(object = output, file = result_file_name)