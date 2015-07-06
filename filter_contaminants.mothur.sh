
# Edit this variable if starting with raw sequence that needs to be dumultiplexed
input_fasta_raw=./test_data/test_seqs_raw.fna

# Edit this variable if starting with fasta in split_libraries.py output format
input_fasta_sl=./test_data/test_seqs.fna


# Output directory
out_dir=./test_output

# Sample metadata mapping filepath
mapping_fp=./test_data/test_seqs_sample_map.txt

# Blank descriptor in sample map ('category:value')
blank_category='Blank:1'

# OTU picking width for clustering Blanks
blank_cluster_pid=.97

# Minimum abundance in Blanks for blank OTU to be considered a potential contaminant
min_blank_pct=.01

# Directory name for filter script output
filter_out_dir=filter_output

# Percent ID for mapping sample sequences against putative blank contaminant OTUs
contaminant_similarity=.97

# Print status messages
verbose=True



mkdir ${out_dir}

#### split libraries ####

# You may need to modify the split libraries command according to your format
if [ ${input_fasta_raw} ];
	then
		if [ ${verbose} = "True" ]; then echo "Splitting libraries"; fi;
		split_libraries.py -m ${mapping_fp} -f ${input_fasta_raw} -o ${out_dir}/split_libraries_out -b 8
		echo "true, splitting libraries"
		input_fasta=${out_dir}/split_libraries_out/seqs.fna
else
		echo "false, setting input fast to slo"
		input_fasta=${input_fasta_sl}
fi



#### Pick unique seqs #####

#pick unique seqs
if [ ${verbose} = "True" ]; then echo "Picking unique seqs"; fi;
pick_otus.py -i ${input_fasta} -s 1 -o ${out_dir}/unique_seqs

# make unique seqs otu table
if [ ${verbose} = "True" ]; then echo "Making unique seq OTU table"; fi;
make_otu_table.py -i ${out_dir}/unique_seqs/*_otus.txt -o ${out_dir}/unique_seqs/unique_seqs_otu_table.biom

# pick unique seqs rep set
if [ ${verbose} = "True" ]; then echo "Picking unique seq rep set"; fi;
pick_rep_set.py -i ${out_dir}/unique_seqs/*_otus.txt -m most_abundant -f ${input_fasta} -o ${out_dir}/unique_seqs/unique_seqs_rep_set.fna


#### Make blank contaminant reference db #####

# filter fasta to just blanks
if [ ${verbose} = "True" ]; then echo "Filtering fasta "; fi;
filter_fasta.py -f ${input_fasta} -o ${out_dir}/blanks.fna --mapping_fp ${mapping_fp} --valid_states ${blank_category}

# cluster blanks
if [ ${verbose} = "True" ]; then echo "Clustering blanks"; fi;
pick_otus.py -i ${out_dir}/blanks.fna -s ${blank_cluster_pid} -o ${out_dir}/blank_clustered

# make blank otu table
if [ ${verbose} = "True" ]; then echo "Making blank OTU table"; fi;
make_otu_table.py -i ${out_dir}/blank_clustered/*_otus.txt -o ${out_dir}/blank_clustered/blanks_otu_table.biom

# pick blank rep set
if [ ${verbose} = "True" ]; then echo "Picking blank rep set"; fi;
pick_rep_set.py -i ${out_dir}/blank_clustered/*_otus.txt -m most_abundant -f ${out_dir}/blanks.fna -o ${out_dir}/blank_clustered/blanks_rep_set.txt

# filter blank otu table
if [ ${verbose} = "True" ]; then echo "Filtering blank OTU table"; fi;
filter_otus_from_otu_table.py -i ${out_dir}/blank_clustered/blanks_otu_table.biom -o ${out_dir}/blank_clustered/blanks_otu_table.filtered.biom --min_count_fraction ${min_blank_pct}

# filter blank rep set
if [ ${verbose} = "True" ]; then echo "Filtering fasta rep set"; fi;
filter_fasta.py -f ${out_dir}/blank_clustered/blanks_rep_set.txt -b ${out_dir}/blank_clustered/blanks_otu_table.filtered.biom  -o ${out_dir}/blank_clustered/blanks_rep_set.filtered.fna

mkdir ${out_dir}/${filter_out_dir}

# run filtering script
if [ ${verbose} = "True" ]; then echo "Running filtering script"; fi;

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
if [ ${verbose} = "True" ]; then echo "Filtering split_library fasta "; fi;
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

if [ ${input_fasta_raw} ]
	then
		# get list of original sequence headers
		if [ ${verbose} = "True" ]; then echo "Getting filtered fasta headers"; fi;
		perl -ne '/^>.+?\s(.+?\n)/ and print $1' ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.fna > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt

		# use seqtk to filter original 
		if [ ${verbose} = "True" ]; then echo "Filtering raw fasta"; fi;
		seqtk subseq ${input_fasta_raw} ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.headers.txt > ${out_dir}/${filter_out_dir}/all_noncontaminant_seqs.raw.fna
fi