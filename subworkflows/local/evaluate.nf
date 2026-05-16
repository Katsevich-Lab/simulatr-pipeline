include { EVALUATE_METHODS } from '../../modules/local/evaluate/main.nf'

workflow EVALUATE {
    take:
    chunk_result
    benchmarking_info

    main:
    chunk_result_grouped     = chunk_result.groupTuple(by: [0, 1])
    benchmarking_info_grouped = benchmarking_info
    evaluate_input            = chunk_result_grouped.join(benchmarking_info_grouped, by: [0, 1])

    EVALUATE_METHODS(evaluate_input)

    emit:
    evaluation_results = EVALUATE_METHODS.out.evaluation_results
    error_json         = EVALUATE_METHODS.out.error_json
    memory_tsv         = EVALUATE_METHODS.out.memory_tsv
}
