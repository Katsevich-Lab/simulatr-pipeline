#!/usr/bin/env Rscript

library(simulatr)

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
simulatr_spec <- readRDS(args[1])
method <- args[2]
row_idx <- as.integer(args[3])
proc_id <- as.integer(args[4])
n_processors <- as.integer(args[5])
B_in <- as.integer(args[6])

# extract data generator and its ordered arguments
data_generator <- simulatr_spec@generate_data_function
ordered_args_data_gen <- get_ordered_args(data_generator, simulatr_spec, row_idx)

# extract the method object and its ordered arguments
method_object <- simulatr_spec@run_method_functions[[method]]
ordered_args_method <- get_ordered_args(method_object, simulatr_spec, row_idx)

if(!data_generator@loop | !method_object@loop){
  # find the row indices
  stop("Data generation or method application with loop = FALSE not yet implemented")
}

# extract the seed
seed <- simulatr_spec@fixed_parameters$seed

# find the replicate indices
B <- if (B_in != 0) B_in else simulatr_spec@fixed_parameters$B
all_b <- 1:B
proc_id_b <- all_b[1 + (all_b %% n_processors) == proc_id]

result_list <- vector(mode = "list", length = length(proc_id_b))
for(b_idx in 1:length(proc_id_b)){
  b <- proc_id_b[b_idx]
  curr_df <- R.utils::withSeed(do.call(data_generator@f, ordered_args_data_gen), seed = seed + b)
  ordered_args_method[[1]] <- curr_df
  out <- dplyr::tibble(output = list(R.utils::withSeed(do.call(method_object@f, ordered_args_method), seed = seed)),
                       run_id = b)
  result_list[[b_idx]] <- out
}
result_df <- do.call(rbind, result_list)

# save result
to_save <- collate_result_list(result_df, proc_id, row_idx, method)
output_filename <- sprintf("chunk_result_%s_%d_%d.rds", 
                           method, row_idx, proc_id)
saveRDS(to_save, output_filename)