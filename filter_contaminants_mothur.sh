# Output directory
out_dir=./test_output/piotr_test_data

mkdir ${out_dir}

# Filter with just Mothur table

python ./qiime_scripts/scripts/decontaminate_script.py \
--mothur_counts_fp ./piotr_script/Test.count_table \
-o ${out_dir} \
--blank_id_fp ./piotr_script/Test_blank_list.txt \
--min_relabund_threshold 0.001 \
--prescreen_threshold 0.7 \
--removal_stat_blank maxB \
--removal_stat_sample maxS \
--removal_differential 10 \
--drop_lib_threshold 0.3 \
--write_per_seq_stats \
--write_per_library_stats \
--write_filtered_output

