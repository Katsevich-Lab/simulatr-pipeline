simspec_filename=$(pwd)/sim_spec_obj.rds
output_dir=$(pwd)
result_filename=test_result.rds
nextflow run main.nf \
    --simulatr_specifier_fp $simspec_filename \
    --result_dir $output_dir \
    --result_file_name $result_filename \
    --B_check 2 \
    --B 20 \
    --max_gb 4 \
    -resume