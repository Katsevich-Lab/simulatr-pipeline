/*
 * SIMULATR workflow.
 *
 * Asserts the installed simulatr package's contract version matches
 * what the pipeline was built for, then fans out:
 *
 *   BENCHMARK -> SIMULATE -> EVALUATE -> COLLECT
 *
 * Memory and error sidecars from each phase land in COLLECT, which
 * publishes memory.parquet / errors.parquet alongside the headline
 * simulatr_result.rds.
 */

include { CHECK_CONTRACT }                  from '../modules/local/check_contract/main.nf'
include { BENCHMARK }                       from '../subworkflows/local/benchmark.nf'
include { SIMULATE }                        from '../subworkflows/local/simulate.nf'
include { EVALUATE }                        from '../subworkflows/local/evaluate.nf'
include { COLLECT }                         from '../subworkflows/local/collect.nf'

workflow SIMULATR {
    CHECK_CONTRACT()

    BENCHMARK()
    SIMULATE(BENCHMARK.out.proc_id_info)
    EVALUATE(SIMULATE.out.chunk_result, BENCHMARK.out.benchmarking_info)

    error_jsons  = BENCHMARK.out.error_json
                    .mix(SIMULATE.out.error_json)
                    .mix(EVALUATE.out.error_json)
    memory_tsvs  = BENCHMARK.out.memory_tsv
                    .mix(SIMULATE.out.memory_tsv)
                    .mix(EVALUATE.out.memory_tsv)

    COLLECT(
        BENCHMARK.out.benchmarking_info,
        EVALUATE.out.evaluation_results,
        error_jsons,
        memory_tsvs
    )
}
