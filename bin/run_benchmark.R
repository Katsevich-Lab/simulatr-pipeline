#!/usr/bin/env Rscript

library(simulatr)

# read command line arguments
args <- commandArgs(trailingOnly = TRUE)
simulatr_spec <- readRDS(args[1])
method <- args[2]
row_idx <- as.integer(args[3])
B_check <- as.integer(args[4])
B_in <- as.integer(args[5])
max_gb <- as.numeric(args[6])
max_hours <- as.numeric(args[7])

# set constants
BYTES_PER_GB <- 2^30
SECONDS_PER_HOUR <- 60 * 60
R_SESSION_GB <- 0.5
HOURS_PROPORTION_USED <- 0.9
GB_PROPORTION_USED <- 0.9

# extract data generator and its ordered arguments
data_generator <- simulatr_spec@generate_data_function
ordered_args_data_gen <- get_ordered_args(data_generator, simulatr_spec, row_idx)

# extract the method object and its ordered arguments
method_object <- simulatr_spec@run_method_functions[[method]]
ordered_args_method <- get_ordered_args(method_object, simulatr_spec, row_idx)

# extract the seed
seed <- simulatr_spec@fixed_parameters$seed

# benchmark data generation
if (data_generator@loop) {
  data_list <- vector(mode = "list", length = B_check)
  data_gb <- numeric(B_check)
  data_hours <- numeric(B_check)
  for(b in 1:B_check){
    benchmark_data <- bench::mark({
      data_list[[b]] <- R.utils::withSeed(
        do.call(data_generator@f, ordered_args_data_gen),
        seed = seed + b
      )
    })
    data_gb[b] <- as.numeric(benchmark_data$mem_alloc) / BYTES_PER_GB
    data_hours[b] <- as.numeric(benchmark_data$median) / SECONDS_PER_HOUR
  }
  data_gb_per_rep <- max(data_gb)
  data_hours_per_rep <- max(data_hours)
} else {
  benchmark_data <- bench::mark(data_list <- do.call(data_generator@f, ordered_args_data_gen))
  data_gb_per_rep <- as.numeric(benchmark_data$mem_alloc) / BYTES_PER_GB / B_check
  data_hours_per_rep <- as.numeric(benchmark_data$median) / SECONDS_PER_HOUR / B_check
}

# benchmark method application
if (method_object@loop){
  result_list <- vector(mode = "list", length = B_check)
  method_gb <- numeric(B_check)
  method_hours <- numeric(B_check)
  for(b in 1:B_check){
    curr_df <- data_list[[b]]
    ordered_args_method[[1]] <- curr_df
    benchmark_method <- bench::mark({
      method_result <- R.utils::withSeed(
        do.call(method_object@f, ordered_args_method),
        seed = seed
      )
    })
    method_result$run_id <- b
    result_list[[b]] <- method_result
    method_gb[b] <- as.numeric(benchmark_method$mem_alloc) / BYTES_PER_GB
    method_hours[b] <- as.numeric(benchmark_method$median) / SECONDS_PER_HOUR
  }
  result_df <- do.call(rbind, result_list)
  method_hours_per_rep <- max(method_hours)
  method_gb_per_rep <- max(method_gb)
} else {
  benchmark_method <- bench::mark(result_df <- do.call(method_object@f, ordered_args_method))
  method_hours_per_rep <- as.numeric(benchmark_method$median) / SECONDS_PER_HOUR / B_check
  method_gb_per_rep <- as.numeric(benchmark_method$mem_alloc) / BYTES_PER_GB / B_check
}
result_gb_per_rep <- as.numeric(pryr::object_size(result_df)) / BYTES_PER_GB / B_check

# compute the number of processors needed
B <- if (B_in != 0) B_in else simulatr_spec@fixed_parameters$B
# maximum number of reps that can be run within the max_gb budget
max_reps_gb <- (max_gb * GB_PROPORTION_USED - R_SESSION_GB - data_gb_per_rep - method_gb_per_rep) / result_gb_per_rep
# maximum number of reps that can be run within the max_hours budget
max_reps_hours <- max_hours * HOURS_PROPORTION_USED / (data_hours_per_rep + method_hours_per_rep)

if(max_reps_gb < 1){
  stop(sprintf("One rep of the %s method requires %f GB of memory, but the 
              maximum memory you specified is %f GB. Please increase the 
              max_gb parameter.", 
              method, 
              R_SESSION_GB + data_gb_per_rep + method_gb_per_rep + result_gb_per_rep, 
              max_gb))
}
if(max_reps_hours < 1){
  stop(sprintf("One rep of the %s method requires %f hours, but the 
              maximum time you specified is %f hours. Please increase the 
              max_hours parameter.", 
              method, 
              data_hours_per_rep + method_hours_per_rep, 
              max_hours))
}

max_reps <- min(max_reps_gb, max_reps_hours)
n_processors <- max(ceiling(B / max_reps))

# write benchmarking information
benchmarking_info <- data.frame(method = method, 
                                grid_id = row_idx, 
                                data_gb_per_rep = data_gb_per_rep,
                                method_gb_per_rep = method_gb_per_rep,
                                result_gb_per_rep = result_gb_per_rep,
                                method_hours_per_rep = method_hours_per_rep, 
                                data_hours_per_rep = data_hours_per_rep,
                                max_reps_gb = max_reps_gb,
                                max_reps_hours = max_reps_hours,
                                max_reps = max_reps,
                                n_processors = n_processors)
benchmarking_info_filename <- sprintf("benchmarking_info_%s_%d.rds", 
                           method, row_idx)
saveRDS(benchmarking_info, benchmarking_info_filename)

# write processors information
proc_id_info <- data.frame(method = method, 
                           grid_id = row_idx, 
                           proc_id = 1:n_processors, 
                           n_processors = n_processors)
proc_id_info_filename <- sprintf("proc_id_info_%s_%d.csv", 
                           method, row_idx)
write.table(proc_id_info, file = proc_id_info_filename, 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ",")