include { OBTAIN_BASIC_INFO; RUN_BENCHMARK } from '../../modules/local/benchmark/main.nf'

workflow BENCHMARK {
    take:
    /* nothing — drives off params */

    main:
    OBTAIN_BASIC_INFO()

    method_names_ch = OBTAIN_BASIC_INFO.out.method_names_raw.splitText().map { it.trim() }
    grid_rows_ch    = OBTAIN_BASIC_INFO.out.grid_rows_raw.splitText().map { it.trim() }
    method_x_grid   = method_names_ch.combine(grid_rows_ch)

    RUN_BENCHMARK(method_x_grid)

    emit:
    proc_id_info       = RUN_BENCHMARK.out.proc_id_info
    benchmarking_info  = RUN_BENCHMARK.out.benchmarking_info
    error_json         = RUN_BENCHMARK.out.error_json
    memory_tsv         = RUN_BENCHMARK.out.memory_tsv
}
