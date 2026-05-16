process OBTAIN_BASIC_INFO {
    tag 'get_info'
    label 'OBTAIN_BASIC_INFO'

    output:
    path 'method_names.txt', emit: method_names_raw
    path 'grid_rows.txt',    emit: grid_rows_raw

    script:
    """
    get_info_for_nextflow.R ${params.simulatr_specifier_fp}
    """
}

process RUN_BENCHMARK {
    tag "method: $method; grid row: $grid_row"
    label 'RUN_BENCHMARK'

    input:
    tuple val(method), val(grid_row)

    output:
    path "proc_id_info_${method}_${grid_row}.csv",                                              emit: proc_id_info
    tuple val(method), val(grid_row), path("benchmarking_info_${method}_${grid_row}.rds"),     emit: benchmarking_info
    path "error_benchmark_${method}_${grid_row}.json",                                          optional: true, emit: error_json
    path '.command.memory.tsv',                                                                 optional: true, emit: memory_tsv

    script:
    """
    bash ${projectDir}/bin/rss_sampler.sh start
    run_benchmark.R \\
        ${params.simulatr_specifier_fp} \\
        ${method} \\
        ${grid_row} \\
        ${params.B_check} \\
        ${params.B} \\
        ${params.max_gb} \\
        ${params.max_hours} \\
        ${params.benchmark_memory} \\
        ${params.R_SESSION_GB} \\
        ${params.GB_PROPORTION_USED}
    bash ${projectDir}/bin/rss_sampler.sh stop
    """
}
