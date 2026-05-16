# Changelog

## v1.0.0 — Phase 3 redesign

* **Restructured under the nf-core directory layout:**
  `workflows/simulatr.nf` top-level workflow, four named subworkflows
  in `subworkflows/local/` (`BENCHMARK`, `SIMULATE`, `EVALUATE`,
  `COLLECT`), and per-process modules in `modules/local/`.
* New typed `nextflow_schema.json` covering `--max_gb`,
  `--max_hours`, `--benchmark_memory`, `--benchmark_memory_base`,
  `--R_SESSION_GB`, `--GB_PROPORTION_USED`, `--array_size`,
  `--simulatr_specifier_fp`, `--result_dir`, `--result_file_name`,
  `--B`, `--B_check`.
* **Slurm array-job support** via `process.array = params.array_size`
  on `RUN_SIMULATION_CHUNK` and `RUN_BENCHMARK`.
* **Layer-2 memory benchmarking.** `process.afterScript` writes
  `sacct`'s `MaxRSS / MaxVMSize` to `.command.memory.tsv` per task on
  Slurm; `bin/rss_sampler.sh` provides the same shape on non-Slurm
  executors. `COLLECT_MEMORY` concatenates into `memory.parquet`.
* **Error-channel JSON sidecars.** Each runner wraps its work in
  `tryCatch` that writes `error_<phase>_<grid_id>_<chunk_id>.json` on
  failure; `COLLECT_ERRORS` concatenates into `errors.parquet`.
* **Pipeline contract assert.** New `CHECK_CONTRACT` process invokes
  `simulatr::contract_version()` at workflow start; mismatches with
  `params.expected_contract_version` (default `1`) abort with a clear
  pointer at the package's `contract.json`.
* **Rewrote `bin/*.R` against the new simulatr verb API.** Phase 2
  removed `simulatr_specifier`, `get_ordered_args`, and friends; the
  runners now use `simulatr_benchmark()`, `simulatr_run()`'s
  interleaved generate-then-run loop, and `simulatr_evaluate()`.
* New `conf/test.config` runs the workflow on `tests/data/study_test.rds`
  (2 grid rows × 1 method × 4 reps) in ~2 minutes on the CI runner.
* `CITATIONS.md` and `docs/cluster.md` added.
* DSL1 fixtures and the legacy spec-builder removed.

## Pre-1.0

* See `git log origin/parallel-evaluation` for the
  `memory-efficient`/`parallel-evaluation` lineage that landed
  decomposed memory accounting and per-(method, grid_row) sharding of
  the evaluation step.
