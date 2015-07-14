
out_dir=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_out
otu_table=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_otu_table.txt
sample_fasta=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_sample_seqs.fna
ref_db=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_contam_ref_db.fna
counts_fp=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_counts_table.txt
mapping_fp=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_sample_map.txt
blanks_fp=./qiime_scripts/qiime_test_data/decontaminate/simple_test/simple_test_blanks.txt
blank_category='Blank:1'

# # Filter with ref and reinstatement

# python ./qiime_scripts/scripts/decontaminate_script.py \
# -i ${otu_table} \
# -o ${out_dir} \
# -m ${mapping_fp} \
# -f ${sample_fasta} \
# -s ${blank_category} \
# --contaminant_db_fp ${ref_db} \
# --contaminant_similarity 0.97 \
# --removal_stat_blank maxB \
# --removal_stat_sample maxS \
# --removal_differential 1 \
# --reinstatement_stat_blank maxB \
# --reinstatement_stat_sample maxS \
# --reinstatement_differential 1 \
# --reinstatement_sample_number 2 \
# --reinstatement_method union 

# Filter with ref and reinstatement

python decontaminate_unitary.py \
--mothur_counts_fp ${counts_fp} \
-o ${out_dir} \
-m ${mapping_fp} \
-f ${sample_fasta} \
-s ${blank_category} \
--max_correlation -0.3 \
--correlate_header Abund \
--write_output_seq_lists \
--min_relabund_threshold 0.2 \
--write_filtered_output \
--drop_lib_threshold 0.8 \
--write_per_library_stats \
--write_per_seq_stats \
--write_per_seq_disposition


# Filter with just Mothur table

# python ./qiime_scripts/scripts/decontaminate_script.py \
# --mothur_counts_fp ${counts_fp} \
# -o ${out_dir} \
# --blank_id_fp ${blanks_fp} \
# --min_relabund_threshold 0.001 \
# --prescreen_threshold 0.5 \
# --removal_stat_blank maxB \
# --removal_stat_sample maxS \
# --removal_differential 1 \




# # filter fasta by all good otu map
# if [ ${verbose} = "True" ]; then echo "Filtering split_library fasta "; fi;
# filter_fasta.py -f ${input_fasta} \
# -o ${out_dir}/${filter_out_dir}/passed_seqs.fna \
# -m ${out_dir}/${filter_out_dir}/passed_otu_map.txt

# # filter fasta by reinstated otu map
# filter_fasta.py -f ${input_fasta} \
# -o ${out_dir}/${filter_out_dir}/reinstated_seqs.fna \
# -m ${out_dir}/${filter_out_dir}/reinstated_contaminants_otu_map.txt

# # cat filtered fastas together into a contaminant-free fasta
# cat ${out_dir}/${filter_out_dir}/passed_seqs.fna ${out_dir}/${filter_out_dir}/reinstated_seqs.fna > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.fna

# ### OPTIONAL STEP: filter original fastq/a by split-libraries-format fasta ####

# if [ ${input_fasta_raw} ]
# 	then
# 		# get list of original sequence headers
# 		if [ ${verbose} = "True" ]; then echo "Getting filtered fasta headers"; fi;
# 		perl -ne '/^>.+?\s(.+?\n)/ and print $1' ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.fna > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt

# 		# use seqtk to filter original 
# 		if [ ${verbose} = "True" ]; then echo "Filtering raw fasta"; fi;
# 		seqtk subseq ${input_fasta_raw} ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.raw.fna
# fi