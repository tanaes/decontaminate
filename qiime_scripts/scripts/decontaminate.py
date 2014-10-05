# Before decontamination:
# remove singletons
# de novo chimera check



# Read in parameters

#!/usr/bin/env python
# File created on 09 Aug 2012

from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jon Sanders", "Nate Bresnick", "Aaron Behr"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"

from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup, make_option
from qiime.parse import parse_qiime_parameters, parse_taxonomy
from qiime.filter import sample_ids_from_metadata_description
from qiime.decontaminate import get_contamination_stats, compare_blank_abundances, print_filtered_otu_map, print_results_file
from biom import load_table
from brokit.usearch import usearch_qf
from brokit.uclust import get_clusters_from_fasta_filepath
import os


options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = """
A script to filter sequences by potential contaminants"""
script_info['script_description'] = """
Description to be added later.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""", """
The following steps are performed by the command below:

1.

2. 
""", """
test_cospeciation_only_nate_draft_2.py -i $PWD/Example_output/cOTUs_test 
-p $PWD/Example_input/otu_table_HostSpecies_rarified_filtered.txt 
-a $PWD/Example_input/host_tree.tre -o $PWD/Example_output/hommola_test 
-T $PWD/hommola -t $PWD/Example_input/taxonomy.txt -m $PWD/Example_input/sample_map.txt 
-c HostSpecies"""))
script_info['output_description'] = """
To be added later
"""
script_info['required_options'] = [
    options_lookup["otu_table_as_primary_input"],
    options_lookup["output_dir"],
    options_lookup["mapping_fp"],
    options_lookup["input_fasta"],
    make_option('-M', '--otu_map_fp', type="existing_filepath",
                 help='the input OTU map file [REQUIRED]'),
    ]
script_info['optional_options'] = [
    make_option('-s',
                '--valid_states', type='string',
                help="string describing valid states (e.g. 'Treatment:Fasting') [default: %default]"),
    make_option('--sample_id_fp',
                type='existing_filepath',
                help='path to file listing sample ids to keep [default: %default]'),
    make_option('--contaminant_db_fp', type="existing_filepath",
              help='A FASTA file of potential contaminant sequences'),
    make_option('-c', '--contaminant_similarity', type='float', default=0.97,
                help=('Sequence similarity threshold for contaminant matches')),
    make_option('--removal_stat_blank', type="choice", choices=["maxB", "avgB"],
                 help='blank statistic to be used for removal (maxB, avgB)'),
    make_option('--removal_stat_sample', type="choice", choices=["maxS", "avgS"],
                 help='sample statistic to be used for removal (maxS, avgS)'),
    make_option('--removal_differential', type="float",
                 help='differential proportion for removal (maxB > X * maxS)'),
    make_option('--reinstatement_stat_blank', type="choice", choices=["maxB", "avgB"],
                 help='blank statistic to be used for reinstatement (maxB, avgB)'),
    make_option('--reinstatement_stat_sample', type="choice", choices=["maxS", "avgS"],
                 help='sample statistic to be used for reinstatement (maxS, avgS)'),
    make_option('--reinstatement_differential', type="float",
                 help='differential proportion for reinstatement (maxS > X * maxB)'),
    make_option('--reinstatement_sample_number', type="int",
                 help='minimum number of samples necessary for reinstatement'),
    make_option('--reinstatement_method', type="choice", choices=["union", "intersection"],
                 help='method to rectify reinstatement criteria')
   


    ]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    valid_states = opts.valid_states
    sample_id_fp = opts.sample_id_fp
    contaminant_db_fp = opts.contaminant_db_fp
    contaminant_similarity = opts.contaminant_similarity
    input_fasta_fp = opts.input_fasta_fp
    otu_map_fp = opts.otu_map_fp
    output_dir = opts.output_dir
    removal_stat_blank = opts.removal_stat_blank
    removal_stat_sample = opts.removal_stat_sample
    removal_differential = opts.removal_differential
    reinstatement_stat_sample = opts.reinstatement_stat_sample
    reinstatement_stat_blank = opts.reinstatement_stat_blank
    reinstatement_differential = opts.reinstatement_differential
    reinstatement_sample_number = opts.reinstatement_sample_number
    reinstatement_method = opts.reinstatement_method

    # Make unique seq OTU table (biom file)

    # Compute unique seq stats
    #   output biom file with unique seq stats

    # Optionally: make candidate contaminant DB
    #   remove sequences present at higher abundance in samples
    #   cluster blanks
    #   remove low-abundance contaminant OTUs

    # Filter by similarity against candidate contaminant DB
    #   annotate unique seq OTU table with top hit (OTU#, rep seq, ID%)
    #   make list of seqs @ threshold

    # Calculate reinstatement rule for filtered sequences

    # Generate lists of seqs failing:
    #   - unique seq rule
    #   - hit to contaminant
    #   - reinstatement after hit

    # open uniqseq biom file
    if not ((mapping_fp and valid_states) or
            sample_id_fp is not None):
        option_parser.error("No filtering requested. Must provide either "
                            "mapping_fp and valid states, min counts, "
                            "max counts, or sample_id_fp (or some combination "
                            "of those).")

    removal_options_counter = 0
    if removal_stat_blank:
        removal_options_counter += 1
    if removal_stat_sample:
        removal_options_counter += 1
    if removal_differential:
        removal_options_counter += 1

    if ((removal_options_counter > 0) and (removal_options_counter < 3)):
        option_parser.error("Must provide all of "
                            "removal_stats_blank, "
                            "removal_stat_sample, and "
                            "removal_differential, or none.")

    reinstatement_options_counter = 0
    if reinstatement_stat_blank:
        reinstatement_options_counter += 1
    if reinstatement_stat_sample:
        reinstatement_options_counter += 1
    if reinstatement_differential:
        reinstatement_options_counter += 1

    if ((reinstatement_options_counter > 0) and (reinstatement_options_counter < 3)):
        option_parser.error("Must provide all of "
                            "reinstatement_stats_blank, "
                            "reinstatement_stat_sample, and "
                            "reinstatement_differential, or none.")

    if (reinstatement_options_counter == 3 and reinstatement_sample_number) and not reinstatement_method:
        option_parser.error("If providing sample number AND abundance criteria "
                            "for sequence reinstatement, must also provide "
                            "a method for combining results.")        

    unique_seq_biom = load_table(otu_table_fp)
        
    # get blank sample IDs from mapping file

    blank_sample_ids = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)

    sample_sample_ids = set(unique_seq_biom.ids(axis='sample')) - set(blank_sample_ids)

    sample_biom = unique_seq_biom.filter(sample_sample_ids, axis='sample', inplace=False)

    # calculate contamination statistics dictionary from uniqseq biom file

    contamination_stats_header, contamination_stats_dict = \
        get_contamination_stats(unique_seq_biom, blank_sample_ids)

    abund_contaminants = set()
    ref_contaminants = set()
    reinstated_seqs = set()

    # Pick seqs that fail the blank vs sample rule
    if removal_stat_blank:
        abund_contaminants = compare_blank_abundances(contamination_stats_dict, 
                                contamination_stats_header,
                                removal_stat_blank,
                                removal_stat_sample,
                                removal_differential,
                                negate=False)

    # Blast against contaminant DB
    if contaminant_db_fp:
        clusters, failures, seeds = get_clusters_from_fasta_filepath(
            input_fasta_fp,
            input_fasta_fp,
            percent_ID=contaminant_similarity,
            max_accepts=1,
            max_rejects=8, 
            stepwords=8,
            word_length=8,
            optimal=False,
            exact=False,
            suppress_sort=False,
            output_dir=output_dir,
            enable_rev_strand_matching=False,
            subject_fasta_filepath=contaminant_db_fp,
            suppress_new_clusters=True,
            return_cluster_maps=True,
            stable_sort=False,
            save_uc_files=True,
            HALT_EXEC=False)

        # Pick seqs that fail the similarity to contaminants rule

        ref_contaminants = set(unique_seq_biom.ids(axis='observation')) - set(failures)

        # Pick seqs from similarity search to reinstate
        if (reinstatement_stat_blank and reinstatement_stat_sample and reinstatement_differential):
            abund_reinstated_seqs = compare_blank_abundances(contamination_stats_dict, 
                                    contamination_stats_header,
                                    reinstatement_stat_sample,
                                    reinstatement_stat_blank,
                                    reinstatement_differential,
                                    negate=False)

            # Only consider seqs as reinstated if previously identified as contaminants
            abund_reinstated_seqs = (ref_contaminants | set(abund_contaminants)) & set(abund_reinstated_seqs)

            if not reinstatement_sample_number:
                reinstated_seqs = abund_reinstated_seqs

        if reinstatement_sample_number:
            incidence_reinstated_seqs = sample_biom.pa().filter(
                    lambda val, id_, metadata: reinstatement_sample_number <= 
                    val.sum(), axis='observation', inplace=False).ids(
                    axis='observation')

            # Only consider seqs as reinstated if previously identified as contaminants
            incidence_reinstated_seqs = (ref_contaminants | set(abund_contaminants)) & set(incidence_reinstated_seqs)

            if not reinstatement_stat_blank:
                reinstated_seqs = incidence_reinstated_seqs

        if reinstatement_sample_number and reinstatement_stat_blank:
            if reinstatement_method == "union":
                reinstated_seqs = abund_reinstated_seqs | incidence_reinstated_seqs
            elif reinstatement_method == "intersection":
                reinstated_seqs = abund_reinstated_seqs & incidence_reinstated_seqs

    # print filtered OTU maps
    ever_good_seqs = set(sample_biom.ids(axis='observation')) - (set(abund_contaminants) | set(ref_contaminants))
    print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'passed_otu_map.txt'), ever_good_seqs)

    if ref_contaminants:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'ref_contaminants_otu_map.txt'), ref_contaminants)

    if abund_contaminants:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'abund_contaminants_otu_map.txt'), abund_contaminants)

    if reinstated_seqs:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'reinstated_contaminants_otu_map.txt'), reinstated_seqs)


    # print log file / per-seq info
    print_results_file(unique_seq_biom.ids(axis='observation'),
                       contamination_stats_header, contamination_stats_dict, 
                       abund_contaminants,
                       ref_contaminants,
                       reinstated_seqs,
                       os.path.join(output_dir,'contamination_summary.txt'))



if __name__ == "__main__":
    main()




