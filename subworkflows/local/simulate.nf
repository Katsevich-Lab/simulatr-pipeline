include { RUN_SIMULATION_CHUNK } from '../../modules/local/simulate/main.nf'

workflow SIMULATE {
    take:
    proc_id_info_csv

    main:
    RUN_SIMULATION_CHUNK(proc_id_info_csv.splitCsv())

    emit:
    chunk_result = RUN_SIMULATION_CHUNK.out.chunk_result
    error_json   = RUN_SIMULATION_CHUNK.out.error_json
    memory_tsv   = RUN_SIMULATION_CHUNK.out.memory_tsv
}
