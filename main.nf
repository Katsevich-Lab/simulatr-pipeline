// Define parameters
params.result_dir = "."
params.result_file_name = "simulatr_result.rds"
params.B = 0
params.B_check = 5
params.max_gb = 8
params.max_hours = 4

// Define processes

process obtain_basic_info {
    tag "get_info"
    memory '2GB'
    time '15m'

    output:
    path "method_names.txt", emit: method_names_raw
    path "grid_rows.txt", emit: grid_rows_raw

    script:
    """
    get_info_for_nextflow.R $params.simulatr_specifier_fp
    """
}

process run_benchmark {
    tag "method: $method; grid row: $grid_row"
    memory '4GB'
    time '2h'

    input:
    tuple val(method), val(grid_row)

    output:
    path "proc_id_info_${method}_${grid_row}.csv", emit: proc_id_info
    tuple val(method), val(grid_row), path("benchmarking_info_${method}_${grid_row}.rds"), emit: benchmarking_info

    script:
    """
    run_benchmark.R $params.simulatr_specifier_fp $method $grid_row $params.B_check $params.B $params.max_gb $params.max_hours
    """
}

process run_simulation_chunk {
    tag "method: $method; grid row: $grid_row; processor: $proc_id"
    memory "$params.max_gb GB"
    time "$params.max_hours h"

    input:
    tuple val(method), val(grid_row), val(proc_id), val(n_processors)

    output:
    tuple val(method), val(grid_row), path("chunk_result_${method}_${grid_row}_${proc_id}.rds"), emit: chunk_result

    script:
    """
    run_simulation_chunk.R $params.simulatr_specifier_fp $method $grid_row $proc_id $n_processors $params.B
    """
}

process evaluate_methods {
    tag "method: $method; grid row: $grid_row"
    maxRetries 6
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    memory { (Math.pow(2, task.attempt - 1) * 6).toInteger() + 'GB' }
    time { (Math.pow(2, task.attempt - 1) * 15).toInteger() + 'm' }

    input:
    tuple val(method), val(grid_row), path('chunk_result*.rds'), path(benchmarking_info)

    output:
    path "${method}_${grid_row}_results.rds", emit: evaluation_results

    script:
    """
    run_evaluation.R $params.simulatr_specifier_fp $method $grid_row chunk_result* benchmarking_info*
    """
}

process collect_results {
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    memory { (Math.pow(2, task.attempt - 1) * 6).toInteger() + 'GB' }
    time { (Math.pow(2, task.attempt - 1) * 15).toInteger() + 'm' }

    input:
    path("*_results.rds")

    publishDir params.result_dir, mode: "copy"

    output:
    path "$params.result_file_name"

    script:
    """
    Rscript -e "
    outputs_list <- list.files(pattern='*_results.rds') |> lapply(readRDS)
    results <- lapply(outputs_list, function(output)(output\\\$results)) |> 
        data.table::rbindlist() |>
        dplyr::as_tibble()
    metrics <- lapply(outputs_list, function(output)(output\\\$metrics)) |> 
        data.table::rbindlist() |>
        dplyr::as_tibble()
    final_output <- list(results = results, metrics = metrics)
    saveRDS(final_output, '${params.result_file_name}')
    "
    """
}

// Workflow definition

workflow {
    obtain_basic_info()
    method_names_ch = obtain_basic_info.out.method_names_raw.splitText().map{it.trim()}
    grid_rows_ch = obtain_basic_info.out.grid_rows_raw.splitText().map{it.trim()}
    method_cross_grid_row_ch = method_names_ch.combine(grid_rows_ch)

    run_benchmark(method_cross_grid_row_ch)
    run_simulation_chunk(run_benchmark.out.proc_id_info.splitCsv())

    chunk_result_grouped = run_simulation_chunk.out.chunk_result.groupTuple(by: [0, 1])
    benchmarking_info_grouped = run_benchmark.out.benchmarking_info
    evaluate_input = chunk_result_grouped.join(benchmarking_info_grouped, by: [0, 1])

    evaluation_results = evaluate_methods(evaluate_input).collect()
    collect_results(evaluation_results)
}