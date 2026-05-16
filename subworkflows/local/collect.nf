include { COLLECT_BENCHMARKING_RESULTS; COLLECT_RESULTS } from '../../modules/local/collect/main.nf'
include { COLLECT_MEMORY }                                from '../../modules/local/collect_memory/main.nf'
include { COLLECT_ERRORS }                                from '../../modules/local/collect_errors/main.nf'

workflow COLLECT {
    take:
    benchmarking_info
    evaluation_results
    error_jsons
    memory_tsvs

    main:
    COLLECT_BENCHMARKING_RESULTS(benchmarking_info.map { it[2] }.collect())
    COLLECT_RESULTS(evaluation_results.collect())
    COLLECT_MEMORY(memory_tsvs.collect(flat: true).ifEmpty([]))
    COLLECT_ERRORS(error_jsons.collect(flat: true).ifEmpty([]))
}
