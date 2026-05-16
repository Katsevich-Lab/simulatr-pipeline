#!/usr/bin/env Rscript
# check_contract.R <expected_schema_version>
#
# Reads inst/extdata/contract.json from the installed simulatr package
# and asserts its schema_version matches <expected_schema_version>. On
# mismatch, prints a clear error pointing at simulatr::contract_version()
# and exits non-zero.

args <- commandArgs(trailingOnly = TRUE)
expected <- as.integer(args[[1L]])

if (!requireNamespace("simulatr", quietly = TRUE)) {
  stop("simulatr is not installed in this R environment. ",
       "Install it from https://github.com/Katsevich-Lab/simulatr ",
       "or use the docker profile.")
}

actual <- simulatr::contract_version()
if (actual != expected) {
  stop(sprintf(
    "Pipeline contract version mismatch: pipeline expected %d, but the installed simulatr package reports %d. ",
    expected, actual),
    "Either upgrade the package (`remotes::install_github('Katsevich-Lab/simulatr@v",
    expected, ".0.0')`) or check out a matching branch of simulatr-pipeline. ",
    "Run `simulatr::contract_version()` in R to inspect the installed version."
  )
}

writeLines(c(
  sprintf("simulatr contract OK: schema_version = %d", actual),
  sprintf("contract_hash = %s", simulatr::contract_hash())
))
