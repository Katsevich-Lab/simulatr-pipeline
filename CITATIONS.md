# `simulatr-pipeline` citations

If you use this pipeline, please cite the underlying R package and the
workflow runner:

## Package

- Katsevich, E., & Barry, T. (2026). *simulatr: Pipeline-First
  Simulation Studies in R.* See `paper.md` in
  [Katsevich-Lab/simulatr](https://github.com/Katsevich-Lab/simulatr)
  for the JoSS submission draft.

## Workflow runner

- Di Tommaso, P., Chatzou, M., Floden, E. W., Prieto Barja, P.,
  Palumbo, E., & Notredame, C. (2017). *Nextflow enables reproducible
  computational workflows.* Nature Biotechnology, 35(4), 316–319.
  doi:[10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

- Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J.,
  Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). *The
  nf-core framework for community-curated bioinformatics pipelines.*
  Nature Biotechnology, 38(3), 276–278.
  doi:[10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)

## Tools used inside the pipeline

The pipeline calls these CRAN packages at runtime; please cite their
upstreams if your work depends on them:

- `bench` — `bench::mark()` for Layer-1 memory accounting.
- `dplyr`, `tibble`, `tidyr`, `purrr` — tidyverse data wrangling.
- `data.table` — fast `rbindlist` for chunk concatenation.
- `R.utils` — `R.utils::withSeed()` for the per-replicate RNG
  protocol.
- `arrow` (optional) — parquet output for `memory.parquet` /
  `errors.parquet`.
