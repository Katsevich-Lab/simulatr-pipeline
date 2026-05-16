process EVALUATE_METHODS {
    tag "method: $method; grid row: $grid_row"
    label 'EVALUATE_METHODS'

    input:
    tuple val(method), val(grid_row), path('chunk_result*.rds'), path(benchmarking_info)

    output:
    path "${method}_${grid_row}_results.rds",                              emit: evaluation_results
    path "error_evaluate_${method}_${grid_row}.json", optional: true,     emit: error_json
    path '.command.memory.tsv',                       optional: true,     emit: memory_tsv

    script:
    """
    bash ${projectDir}/bin/rss_sampler.sh start
    run_evaluation.R \\
        ${params.simulatr_specifier_fp} \\
        ${method} \\
        ${grid_row} \\
        chunk_result* \\
        ${benchmarking_info}
    bash ${projectDir}/bin/rss_sampler.sh stop
    """
}
