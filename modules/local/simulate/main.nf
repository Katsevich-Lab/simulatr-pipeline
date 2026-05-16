process RUN_SIMULATION_CHUNK {
    tag "method: $method; grid row: $grid_row; processor: $proc_id"
    label 'RUN_SIMULATION_CHUNK'

    input:
    tuple val(method), val(grid_row), val(proc_id), val(n_processors)

    output:
    tuple val(method), val(grid_row), path("chunk_result_${method}_${grid_row}_${proc_id}.rds"), emit: chunk_result
    path "error_chunk_${method}_${grid_row}_${proc_id}.json",                                   optional: true, emit: error_json
    path '.command.memory.tsv',                                                                 optional: true, emit: memory_tsv

    script:
    """
    bash ${projectDir}/bin/rss_sampler.sh start
    run_simulation_chunk.R \\
        ${params.simulatr_specifier_fp} \\
        ${method} \\
        ${grid_row} \\
        ${proc_id} \\
        ${n_processors} \\
        ${params.B}
    bash ${projectDir}/bin/rss_sampler.sh stop
    """
}
