#!/usr/bin/env Rscript
# Builds the small simulatr_study used by `nextflow run . -profile test`.
# Writes tests/data/study_test.rds.
#
# Run with:  Rscript tests/build_test_fixture.R

suppressPackageStartupMessages(library(simulatr))

parameter_grid <- data.frame(
  n        = c(20, 40),
  p        = 5,
  s        = 2,
  beta_val = 1
)

get_ground_truth <- function(p, s, beta_val) {
  beta <- numeric(p)
  beta[seq_len(s)] <- beta_val
  list(beta = beta)
}
parameter_grid <- add_ground_truth(parameter_grid, get_ground_truth)

lm_model <- simulatr_model(
  name     = "lm_model",
  generate = function(n, p, ground_truth) {
    X <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(X %*% ground_truth$beta + rnorm(n))
    list(X = X, y = y)
  }
)

ols <- simulatr_method(
  name = "ols",
  run  = function(data) list(beta = unname(lm.fit(data$X, data$y)$coefficients))
)

rmse <- simulatr_metric(
  name  = "rmse",
  score = function(output, ground_truth) sqrt(mean((output$beta - ground_truth$beta) ^ 2))
)

study <- simulatr_study(
  parameter_grid   = parameter_grid,
  fixed_parameters = list(B = 4, seed = 1),
  models  = list(lm_model = lm_model),
  methods = list(ols = ols),
  metrics = list(rmse = rmse)
)

dir.create("tests/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(study, "tests/data/study_test.rds")
cat("Wrote tests/data/study_test.rds\n")
