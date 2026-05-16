#!/usr/bin/env Rscript
# get_info_for_nextflow.R <study.rds>
#
# Inspect a simulatr_study and emit the cartesian product of method
# names and grid row indices to two text files (one value per line),
# which Nextflow turns into channels.

library(simulatr)

args <- commandArgs(trailingOnly = TRUE)
study <- readRDS(args[[1L]])

method_names <- names(study@methods)
grid_ids     <- seq_len(nrow(study@parameter_grid))

write_vector <- function(file_name, vec) {
  con <- file(file_name)
  writeLines(as.character(vec), con)
  close(con)
}

write_vector("method_names.txt", method_names)
write_vector("grid_rows.txt",    grid_ids)
