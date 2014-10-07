
input_fasta_raw=./test_data/test_seqs_raw.fna
out_dir=test_output
mapping_fp=./test_data/test_seqs_sample_map.txt
blank_category='Blank:1'
blank_cluster_pid=.97
min_blank_pct=.01
filter_out_dir=filter_output
contaminant_similarity=.97

mkdir ${out_dir}


#### split libraries ####
split_libraries.py -m ${mapping_fp} -f ${input_fasta_raw} -o ${test_output}/split_libraries_out -b 8

input_fasta=${test_output}/split_libraries_out/seqs.fna

#### Pick unique seqs #####

#pick unique seqs
pick_otus.py -i ${input_fasta} -s 1 -o ${out_dir}/unique_seqs

# make unique seqs otu table
make_otu_table.py -i ${out_dir}/unique_seqs/*_otus.txt -o ${out_dir}/unique_seqs/unique_seqs_otu_table.biom

# pick unique seqs rep set
pick_rep_set.py -i ${out_dir}/unique_seqs/*_otus.txt -m most_abundant -f ${input_fasta} -o ${out_dir}/unique_seqs/unique_seqs_rep_set.fna


#### Make blank contaminant reference db #####

# filter fasta to just blanks
filter_fasta.py -f ${input_fasta} -o ${out_dir}/blanks.fna --mapping_fp ${mapping_fp} --valid_states ${blank_category}

# cluster blanks
pick_otus.py -i ${out_dir}/blanks.fna -s ${blank_cluster_pid} -o ${out_dir}/blank_clustered

# make blank otu table
make_otu_table.py -i ${out_dir}/blank_clustered/*_otus.txt -o ${out_dir}/blank_clustered/blanks_otu_table.biom

# pick blank rep set
pick_rep_set.py -i ${out_dir}/blank_clustered/*_otus.txt -m most_abundant -f ${out_dir}/blanks.fna -o ${out_dir}/blank_clustered/blanks_rep_set.txt

# filter blank otu table
filter_otus_from_otu_table.py -i ${out_dir}/blank_clustered/blanks_otu_table.biom -o ${out_dir}/blank_clustered/blanks_otu_table.filtered.biom --min_count_fraction ${min_blank_pct}

# filter blank rep set
filter_fasta.py -f ${out_dir}/blank_clustered/blanks_rep_set.txt -b ${out_dir}/blank_clustered/blanks_otu_table.filtered.biom  -o ${out_dir}/blank_clustered/blanks_rep_set.filtered.fna

mkdir ${out_dir}/${filter_out_dir}

# run filtering script
python ./decontaminate.py \
-i ${out_dir}/unique_seqs/unique_seqs_otu_table.biom \
-o ${out_dir}/${filter_out_dir} \
-m ${mapping_fp} \
-f ${out_dir}/unique_seqs/unique_seqs_rep_set.fna \
-M ${out_dir}/unique_seqs/*_otus.txt \
-s ${blank_category} \
--contaminant_db_fp ${out_dir}/blank_clustered/blanks_rep_set.filtered.fna \
--contaminant_similarity ${contaminant_similarity} \
--removal_stat_blank maxB \
--removal_stat_sample maxS \
--removal_differential 1 \
--reinstatement_stat_blank maxB \
--reinstatement_stat_sample maxS \
--reinstatement_differential 1 \
--reinstatement_sample_number 2 \
--reinstatement_method intersection 

# filter fasta by all good otu map
filter_fasta.py -f ${input_fasta} \
-o ${out_dir}/${filter_out_dir}/passed_seqs.fna \
-m ${out_dir}/${filter_out_dir}/passed_otu_map.txt

# filter fasta by reinstated otu map
filter_fasta.py -f ${input_fasta} \
-o ${out_dir}/${filter_out_dir}/reinstated_seqs.fna \
-m ${out_dir}/${filter_out_dir}/reinstated_contaminants_otu_map.txt

# cat filtered fastas together into a contaminant-free fasta
cat ${out_dir}/${filter_out_dir}/passed_seqs.fna ${out_dir}/${filter_out_dir}/reinstated_seqs.fna > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.fna

### OPTIONAL STEP: filter original fastq/a by split-libraries-format fasta ####

# get list of original sequence headers
perl -ne '/^>.+?\s(.+?\n)/ and print $1' ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.fna > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt

# use seqtk to filter original 
seqtk subseq ${input_fasta_raw} ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.raw.fna
