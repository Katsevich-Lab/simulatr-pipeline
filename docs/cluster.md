# Running `simulatr-pipeline` on a cluster

This page covers Slurm specifically; LSF / SGE work analogously after
swapping `process.executor = 'lsf'` / `'sge'` in `conf/slurm.config`.

## Quick start

```bash
nextflow run Katsevich-Lab/simulatr-pipeline \
    -profile slurm \
    --simulatr_specifier_fp /path/to/study.rds \
    --max_gb 16 \
    --max_hours 4 \
    --array_size 200
```

## What the Slurm profile does

`conf/slurm.config` does two things on top of `conf/base.config`:

### 1. Array submission

```nextflow
process {
    withName: 'RUN_SIMULATION_CHUNK' { array = params.array_size }
    withName: 'RUN_BENCHMARK'        { array = params.array_size }
}
```

Without arrays, a 20-method × 4-row × 8-processor study issues 640
`sbatch` calls. With `array_size = 200`, the same study collapses
into 4 array submissions. Scheduler overhead drops accordingly.

The default `200` matches most Slurm clusters' `MaxArraySize`; lower
it if your cluster's limit is smaller.

### 2. Layer-2 memory accounting

```nextflow
process {
    afterScript = '''
        if [ -n "${SLURM_JOB_ID:-}" ]; then
            sacct -j "${SLURM_JOB_ID}" \
                  --format JobID,JobName,MaxRSS,MaxVMSize,Elapsed \
                  --noheader --parsable2 > .command.memory.tsv 2>/dev/null
        fi
    '''
}
```

`sacct` is the cluster's bookkeeping. Per-process peak RSS lands in
`.command.memory.tsv` per task; the `COLLECT_MEMORY` process
concatenates these into `memory.parquet`. The Layer-1 numbers from
`bench::mark()` are already in `benchmarking_results.rds` — diff the
two to find where R's view of memory under-counts (typically Rcpp
arenas, BLAS scratch, or anything called via `reticulate`).

## Sizing resources

The `--max_gb` and `--max_hours` flags constrain *one* task. The
benchmarking phase decides how many `n_processors` per
`(method, grid_row)` fit inside that budget by

```
n_processors = ceil(B / min(max_reps_gb, max_reps_hours))
```

where `max_reps_gb = (max_gb * GB_PROPORTION_USED - R_SESSION_GB
- data_gb_per_rep - method_gb_per_rep) / result_gb_per_rep`.

If `max_reps_gb < 1`, the pipeline pre-stops with a message naming
the offending `(method, grid_row)` — increase `--max_gb` or split the
parameter grid.

## Retries

`conf/base.config` sets `errorStrategy = { task.attempt < 6 && ... }`
with exponential resource doubling (`memory = base * 2^(attempt-1)`).
A first attempt with `--max_gb 16` retries at 32 GB, 64 GB, ..., up
to attempt 6. The same doubling applies to `time`.

## Disabling memory benchmarking

On hosts where `bench::mark` mis-reports (some Docker images), set
`--benchmark_memory 0`. This treats `max_gb` as `Inf` for the purposes
of `n_processors` computation. Layer-2 (`sacct`) numbers are still
collected.
