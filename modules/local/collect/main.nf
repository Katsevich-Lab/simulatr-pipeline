process COLLECT_BENCHMARKING_RESULTS {
    label 'COLLECT_BENCHMARKING_RESULTS'

    input:
    path benchmarking_info

    publishDir params.result_dir, mode: 'copy'

    output:
    path 'benchmarking_results.rds'

    script:
    """
    Rscript -e '
    list.files() |>
      lapply(readRDS) |>
      data.table::rbindlist() |>
      dplyr::as_tibble() |>
      saveRDS(file = "benchmarking_results.rds")
    '
    """
}

process COLLECT_RESULTS {
    label 'COLLECT_RESULTS'

    input:
    path '*_results.rds'

    publishDir params.result_dir, mode: 'copy'

    output:
    path "${params.result_file_name}"

    script:
    """
    Rscript -e "
    outputs_list <- list.files(pattern='*_results.rds') |> lapply(readRDS)
    results <- lapply(outputs_list, function(o) o\\\$results) |>
      data.table::rbindlist() |> dplyr::as_tibble()
    metrics <- lapply(outputs_list, function(o) o\\\$metrics) |>
      data.table::rbindlist() |> dplyr::as_tibble()
    out <- list(
      results = results,
      metrics = metrics,
      session_info = simulatr::simulatr_session_info()
    )
    saveRDS(out, '${params.result_file_name}')
    "
    """
}
