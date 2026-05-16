# `simulatr-pipeline`

[![CI](https://github.com/Katsevich-Lab/simulatr-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/Katsevich-Lab/simulatr-pipeline/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-template-1A9655.svg)](https://nf-co.re/)

Nextflow workflow that runs [`simulatr`](https://github.com/Katsevich-Lab/simulatr)
studies at scale — laptop, Docker, or HPC. The pipeline benchmarks
per-grid-row memory and runtime, shards work into Slurm array jobs,
retries on out-of-memory, and emits memory/error sidecars alongside
the headline result.

## Quick start

```bash
# Build a small study fixture (R + the simulatr package required)
Rscript tests/build_test_fixture.R

# Run end-to-end with the test profile (~2 min on a laptop)
nextflow run . -profile test
```

The published outputs land in `results-test/`:

```
results-test/
├── simulatr_result.rds         # combined results + metrics
├── benchmarking_results.rds    # per-(method, grid_row) Layer-1 bench
├── memory.parquet              # Layer-2 (sacct or rss_sampler) memory
└── errors.parquet              # one row per failed task (if any)
```

## Layout (nf-core)

```
.
├── main.nf                     # entrypoint — calls workflows/simulatr.nf
├── workflows/simulatr.nf       # SIMULATR top-level workflow
├── subworkflows/local/         # BENCHMARK / SIMULATE / EVALUATE / COLLECT
├── modules/local/<name>/       # one process per module
├── conf/
│   ├── base.config             # per-process resource defaults
│   ├── test.config             # minimal end-to-end fixture
│   ├── docker.config           # container = ghcr.io/katsevich-lab/simulatr
│   └── slurm.config            # array submission + sacct memory hook
├── nextflow_schema.json        # typed parameter schema
├── bin/                        # R + bash runners (executed by modules)
└── tests/build_test_fixture.R  # builds tests/data/study_test.rds
```

## Running on Slurm

```bash
nextflow run Katsevich-Lab/simulatr-pipeline \
    -profile slurm \
    --simulatr_specifier_fp /path/to/study.rds \
    --max_gb 16 \
    --max_hours 4 \
    --array_size 200
```

The `slurm` profile (`conf/slurm.config`) does two things:

1. **Array submission.** `process.array = params.array_size` collapses
   many `RUN_SIMULATION_CHUNK` and `RUN_BENCHMARK` tasks into a single
   Slurm array job per process. Reduces scheduler overhead from one job
   per `(method, grid_row, chunk)` to one array per process.
2. **Layer-2 memory accounting.** `process.afterScript` runs
   `sacct -j $SLURM_JOB_ID --format MaxRSS,MaxVMSize` and writes
   `.command.memory.tsv` per task. The `COLLECT_MEMORY` process
   concatenates them into `memory.parquet`.

See `docs/cluster.md` for further detail.

## Package ↔ pipeline contract

The package ships `inst/extdata/contract.json` describing the shard
layout, parquet schemas, and RNG protocol. At workflow start,
`CHECK_CONTRACT` calls `simulatr::contract_version()` and asserts it
matches `params.expected_contract_version` (default `1`). On mismatch,
the workflow exits non-zero with a clear error pointing the user at
how to align the two.

## Errors

Each `bin/*.R` runner wraps its work in `tryCatch` that on error
writes `error_<phase>_<grid_id>_<chunk_id>.json` to the task's work
directory with `{phase, grid_id, chunk_id, run_id, type, message,
traceback}`, then re-throws. The `COLLECT_ERRORS` process concatenates
those JSONs into `errors.parquet`. No more `grep`-ing through
`.command.err` files.

## Memory accounting

Two layers land in `memory.parquet`:

- **Layer 1 (R-heap):** `simulatr_benchmark()` runs `bench::mark`
  during the benchmarking phase; numbers go into
  `benchmarking_results.rds`.
- **Layer 2 (kernel RSS):** Slurm `sacct` (or the portable
  `bin/rss_sampler.sh` for non-Slurm executors) dumps
  `MaxRSS / MaxVMSize` per task into `.command.memory.tsv`, which
  `COLLECT_MEMORY` rolls up into `memory.parquet`.

See the package vignette
`vignette("reading-your-memory-profile", package = "simulatr")` for how
to diff the two.

## See also

- [`simulatr`](https://github.com/Katsevich-Lab/simulatr) (R package).
- [`paper.md`](https://github.com/Katsevich-Lab/simulatr/blob/main/paper.md) (JoSS submission draft).
- [Nextflow](https://www.nextflow.io/) and the [nf-core
  community](https://nf-co.re/).
