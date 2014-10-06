#!/usr/bin/env python
# File created on 09 Aug 2012

from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"

from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup, make_option
from qiime.parse import parse_qiime_parameters, parse_taxonomy
from qiime.filter import sample_ids_from_metadata_description
#from qiime.decontaminate import get_contamination_stats, compare_blank_abundances, print_filtered_otu_map, print_results_file

from biom.parse import parse_biom_table

from qiime.pycogent_backports.uclust import get_clusters_from_fasta_filepath
from qiime.pycogent_backports.usearch import usearch_qf

import os
import numpy as np


options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = """
A script to filter sequences by potential contaminants"""
script_info['script_description'] = """
This script performs a series of filtering steps on a sequence file with the 
intent of removing contaminant sequences. It requires input of an OTU table, a
sample map, an OTU map, a sequence FASTA file, and an output directory.

There are two primary approaches the script can take: (1) comparing sequence
abundances in blank control sequence libraries to those in sample libraries,
where sequences present in blanks are presumed to be contaminants, and (2)
comparing sequences in sample libraries to a database of known contaminants.

In approach (1), OTUs (or unique sequences, if OTU table and map are defined at
100% identity) are tested for their maximum and mean presence in blank and 
sample libraries, and excluded if they satisfy the given criteria. For example, 
if you want to exclude any sequences whose maximum abundance in a blank sample
is more than 10% the maximum abundance in a sample (maxB > 0.1 * maxS), you
would choose '--removal_stat_blank maxB --removal_stat_sample maxS 
--removal_differential 0.1'. For this approach, you must also provide a column 
in your mapping file that indicates which samples to use as blanks, and pass
this information to the script with the 'valid states' option (e.g. 
'Blank:True')

In approach (2), you must provide a fasta library of putative contaminants.
These may be previously clustered OTUs from the blank samples, commonly
sequenced contaminants (if known), or another fasta file. Sequences will be
clustered against this fasta file using Uclust-Ref, and any that match within
a given percent similarity (using the '-c' or '--contaminant_similarity' option)
will be marked as putative contaminants. 

When using approach (2), it is possible to remove 'real' sequences from samples
that just happen to be similar to contaminants. This may be detectable when
using unique sequence OTU tables/maps as input, if the 'real' sequences are
nonetheless slightly different from contaminants. In this case, it may be
desireable to reinstate those unique sequences that are present in samples but
not in blanks. You may do this using criteria of relative abundance (similar to
approach [1], where a sequence is reinstated if its max presence in a sample is
greater than its max presence in a blank, i.e. maxS > X * maxB) or of incidence
in non-blank samples (i.e. reinstated if present in two or more samples). If
both criteria are provided, you must choose to reinstate either the intersection
of the criteria (i.e. BOTH more abundant in samples AND present in 2 or more)
or the union (i.e. EITHER more abundant in samples OR present in 2 or more).
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""", """
The following steps are performed by the command below:

1. Calculate max relative abundance of each sequence in samples and blanks

2. Identify sequences whose maximum abunance in blanks is more than 10% their
maximum abundance in samples.

3. Output OTU maps of sequences for which above is true, and for which above is
false.
""", """
decontaminate.py -i unique_seqs_otu_table.biom -o filter_out_dir 
-m metadata_mapping_file.txt -f unique_seqs_rep_set.fna 
-M unique_seqs_otus.txt -s 'Blank:True' --removal_stat_blank maxB 
--removal_stat_sample maxS --removal_differential 0.1
"""))
script_info['output_description'] = """
This script will output a tab-delimited summary table, indicating the relative
abundance stats for each sequence considered, along with its fate at each step
of the process. 

It will also output an OTU map for each category of sequences identified (e.g.
those never identified as contaminants, those identified as reference-based
contaminants, those identified as abundance-based contaminants, and those
reinstated). These OTU maps can then be used to filter in the input FASTA file. 

Output file naming:
contamination_summary.txt -- tab-delimited per-sequence summary file
assed_otu_map.txt -- OTU map of non-contaminant sequences
ref_contaminants_otu_map.txt -- OTU map of reference contaminant sequences
abund_contaminants_otu_map.txt -- OTU map of abundance contaminant sequences
reinstated_contaminants_otu_map.txt -- OTU map of reinstated sequences
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


def get_contamination_stats(biom_file, blank_sample_ids, proportional=False):
    if not proportional:
        biom_file = biom_file.normObservationBySample()

    blank_data = biom_file.filterSamples(lambda val, id_, metadata: 
        id_ in blank_sample_ids, invert=False)

    sample_data = biom_file.filterSamples(lambda val, id_, metadata: 
        id_ in blank_sample_ids, invert=True)

    stats_dict = {}

    for otu in biom_file.ObservationIds:

        stats_dict[otu] = [sample_data.observationData(otu).max(),
                           sample_data.observationData(otu).mean(),
                           blank_data.observationData(otu).max(),
                           blank_data.observationData(otu).mean()]

    return(['maxS','avgS','maxB','avgB'], stats_dict)

def compare_blank_abundances(stats_dict, stats_header,
                            sample_stat, blank_stat, scalar=1, negate=False):
    
    i_s = stats_header.index(sample_stat)
    i_b = stats_header.index(blank_stat)

    passed_otus = set()

    for otu in stats_dict:
        if(float(stats_dict[otu][i_s]) > (float(scalar) * float(stats_dict[otu][i_b]))):
            passed_otus.add(otu)

    # print passed_otus
    return(passed_otus)

def print_filtered_otu_map(input_otu_map_fp, output_otu_map_fp, filter_set):

    results = set()
    output_otu_map_f = open(output_otu_map_fp, 'w')
    for line in open(input_otu_map_fp, 'U'):
        seq_identifier = line.strip().split('\t')[0]
        # only write this line if the otu has more than n sequences (so
        # greater than n tab-separated fields including the otu identifier)
        if seq_identifier in filter_set:
            output_otu_map_f.write(line)
    output_otu_map_f.close()
    return

def print_results_file(seq_list, 
                       stats_header, stats_dict, 
                       abund_contaminants,
                       ref_contaminants,
                       reinstated_seqs, output_fp):
    
    output_f = open(output_fp, 'w')
    header = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format('SeqID',
             stats_header[0],stats_header[1],stats_header[2],stats_header[3],
             'Abund_contaminant','Ref_contaminant','Reinstated')
    output_f.write(header)
    for otu in seq_list:
        outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(str(otu),
             stats_dict[otu][0],stats_dict[otu][1],stats_dict[otu][2],stats_dict[otu][3],
             1 if otu in abund_contaminants else 0,
             1 if otu in ref_contaminants else 0,
             1 if otu in reinstated_seqs else 0)
        output_f.write(outline)

    return


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

    unique_seq_biom = parse_biom_table(open(otu_table_fp,'Ur'))
        
    # get blank sample IDs from mapping file

    blank_sample_ids = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)

    sample_sample_ids = set(unique_seq_biom.SampleIds) - set(blank_sample_ids)

    sample_biom = unique_seq_biom.filterSamples(lambda val, id_, metadata: 
        id_ in blank_sample_ids, invert=True)

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

        ref_contaminants = set(unique_seq_biom.ObservationIds) - set(failures)

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
            incidence_reinstated_seqs_idx = sample_biom.nonzeroCounts(
                'observation', binary=True) >= reinstatement_sample_number
            incidence_reinstated_seqs = np.array(sample_biom.ObservationIds)[
                                                 incidence_reinstated_seqs_idx]

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
    ever_good_seqs = set(sample_biom.ObservationIds) - (set(abund_contaminants) | set(ref_contaminants))
    print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'passed_otu_map.txt'), ever_good_seqs)

    if ref_contaminants:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'ref_contaminants_otu_map.txt'), ref_contaminants)

    if abund_contaminants:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'abund_contaminants_otu_map.txt'), abund_contaminants)

    if reinstated_seqs:
        print_filtered_otu_map(otu_map_fp, os.path.join(output_dir,'reinstated_contaminants_otu_map.txt'), reinstated_seqs)


    # print log file / per-seq info
    print_results_file(unique_seq_biom.ObservationIds,
                       contamination_stats_header, contamination_stats_dict, 
                       abund_contaminants,
                       ref_contaminants,
                       reinstated_seqs,
                       os.path.join(output_dir,'contamination_summary.txt'))



if __name__ == "__main__":
    main()




