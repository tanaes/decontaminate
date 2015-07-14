#!/usr/bin/env python
# File created on 09 Aug 2012

from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"


from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup, make_option
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_mapping_file_to_dict
from qiime.filter import sample_ids_from_metadata_description
from bfillings.uclust import get_clusters_from_fasta_filepath
from bfillings.usearch import usearch_qf

from scipy.stats import spearmanr
import os.path

from biom import load_table

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
    options_lookup["output_dir"]
    ]
script_info['optional_options'] = [
    options_lookup["otu_table_as_primary_input"],
    make_option('--mothur_counts_fp',
                type='existing_filepath',
                help='path to mothur counts table as input'),
    options_lookup["mapping_fp"],
    make_option('-M', '--otu_map_fp', type="existing_filepath",
                 help='the input OTU map file'),
    make_option('-s',
                '--valid_states', type='string',
                help="Column header:value pair in mapping file identifying blank samples"),
    make_option('--blank_id_fp',
                type='existing_filepath',
                help='path to file listing blank sample ids'),
    options_lookup["input_fasta"],
    make_option('--contaminant_db_fp', type="existing_filepath",
              help='A FASTA file of potential contaminant sequences'),
    make_option('-c', '--contaminant_similarity', type='float', default=0.97,
                help=('Sequence similarity threshold for contaminant matches')),
    make_option('-r', '--max_correlation', type='float',
                help=('Maximum Spearman correlation for contaminant identification')),
    make_option('--correlate_header', type='string',
                help=('Column header in mapping file with correlation data')),
    make_option('--min_relabund_threshold', type="float",
                help='discard sequences below this relative abundance threshold'),
    make_option('--prescreen_threshold', type="float",
                help='prescreen libraries that lose more than this proportion of sequences'),
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
                 help='method to rectify reinstatement criteria'),
    make_option('--drop_lib_threshold', type="float", 
                 help='read loss threshold to drop libraries from output table'),
    make_option('--write_filtered_output', action="store_true", 
                 help='write an output table filtered of contaminants'),
    make_option('--write_per_library_stats', action="store_true", 
                 help='write a per-library decontamination summary'),
    make_option('--write_per_seq_stats', action="store_true", 
                 help='write a per-sequence decontamination summary'),
    make_option('--write_per_seq_disposition', action="store_true", 
                 help='write a per-sequence disposition file'),
    make_option('--write_output_seq_lists', action="store_true",
                 help='write separate sequence name lists for each contaminant category')

    ]

script_info['version'] = __version__


def pick_ref_contaminants(queries, ref_db_fp, input_fasta_fp, contaminant_similarity, output_dir):
    # Blast against contaminant DB

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
        subject_fasta_filepath=ref_db_fp,
        suppress_new_clusters=True,
        return_cluster_maps=True,
        stable_sort=False,
        save_uc_files=True,
        HALT_EXEC=False)

    # Pick seqs that fail the similarity to contaminants rule

    ref_contaminants = set(queries) - set(failures)

    return(ref_contaminants)


def pick_corr_contaminants(sample_biom,
                           corr_data_dict,
                           max_r):
    
    # Filter biom to only samples for which correlate data available
    sample_biom_filt = sample_biom.filter(
        lambda val, id_, metadata: id_ in corr_data_dict, 
        invert=False,
        inplace=False)

    otus = sample_biom_filt.ids(axis='observation')
    samples = sample_biom_filt.ids(axis='sample')

    # Make array of correlate data in same order as biom file
    correlate = [corr_data_dict[x] for x in samples]

    obs_corr_dict = {}

    # Make a 2D array of normalized biom table values
    norm_array = sample_biom_filt.norm(inplace=False).matrix_data.toarray()

    t = 0

    for otu in otus:
        obs_corr_dict[otu] = spearmanr(norm_array[t], correlate)

        t += 1

    # get keys (otu names) for OTUs with less than minimum correlation
    obs_corr_contaminants = [x for x in obs_corr_dict if obs_corr_dict[x][0] < max_r]

    return(set(obs_corr_contaminants), obs_corr_dict)

def reinstate_abund_seqs(putative_contaminants, 
                         contamination_stats_dict, 
                         contamination_stats_header,
                         reinstatement_stat_sample,
                         reinstatement_stat_blank,
                         reinstatement_differential):

    abund_reinstated_seqs = compare_blank_abundances(contamination_stats_dict, 
                            contamination_stats_header,
                            reinstatement_stat_sample,
                            reinstatement_stat_blank,
                            reinstatement_differential,
                            negate=False)
    
    # Only consider seqs as reinstated if previously identified as contaminants
    abund_reinstated_seqs = set(putative_contaminants) & set(abund_reinstated_seqs)

    return(abund_reinstated_seqs)

def reinstate_incidence_seqs(putative_contaminants,
                         unique_seq_biom,
                         blank_sample_ids,
                         reinstatement_sample_number):
    
    sample_biom = unique_seq_biom.filter(lambda val, id_, metadata: 
        id_ in blank_sample_ids, invert=True, inplace=False)

    incidence_reinstated_seqs = sample_biom.pa().filter(
            lambda val, id_, metadata: val.sum() >= reinstatement_sample_number,
            axis='observation', inplace=False).ids(
            axis='observation')

    # Only consider seqs as reinstated if previously identified as contaminants
    incidence_reinstated_seqs = set(putative_contaminants) & set(incidence_reinstated_seqs)

    return(incidence_reinstated_seqs)

def mothur_counts_to_biom(mothur_fp):

    mothur_biom = load_table(mothur_fp)
    mothur_biom.type = u'OTU table'
    filter_biom = mothur_biom.filter(
        lambda val, id_, metadata: id_ in 'total', invert=True)

    return(filter_biom)


def biom_to_mothur_counts(biom_obj):

    sample_ids = biom_obj.ids(axis='sample')
    otu_ids = biom_obj.ids(axis='observation')
    otu_totals = biom_obj.sum(axis='observation')

    outstring = 'Representative_Sequence\ttotal\t' + '\t'.join(sample_ids) + '\n'

    for otu in otu_ids:
        otu_data = biom_obj.data(id = otu, axis = 'observation')
        outstring += '{0}\t{1}\t{2}\n'.format(otu,
                                              int(otu_data.sum()),
                                              '\t'.join(str(x) for x in otu_data.astype('int')))

    return(outstring)




def prescreen_libraries(unique_seq_biom,
                        blank_sample_ids,
                        removal_stat_sample, 
                        removal_stat_blank, 
                        removal_differential, 
                        prescreen_threshold):

    contamination_stats_header, contamination_stats_dict = \
            get_contamination_stats(unique_seq_biom, blank_sample_ids)

    abund_contaminants = compare_blank_abundances(contamination_stats_dict, 
                                contamination_stats_header,
                                removal_stat_sample,
                                removal_stat_blank,
                                removal_differential,
                                negate=True)

    # make relabund table
    norm_biom = unique_seq_biom.norm(inplace = False)

    # filter out sequences marked as contaminants
    norm_biom.filter(lambda val, id_, metadata: id_ in abund_contaminants,
                     axis='observation', invert=True, inplace=True)

    # filter out samples above threshold
    norm_biom.filter(lambda val, id_, metadata: sum(val) > prescreen_threshold,
                          axis='sample', invert=False, inplace=True)

    # Now only have samples failing the prescreening
    above_threshold_samples = norm_biom.ids(axis='sample')

    return above_threshold_samples


def get_contamination_stats(biom_file, blank_sample_ids=None, exp_sample_ids=[], proportional=False):
    if not proportional:
        biom_file = biom_file.norm(inplace=False)

    header = ['maxS','avgS']

    # Calculate blank stats if blank sample names are provided
    if blank_sample_ids:
        blanks = True

        blank_data = biom_file.filter(blank_sample_ids, axis='sample', 
                                      invert=False, inplace=False).matrix_data
        maxB = blank_data.max(axis=1).todense().tolist()
        avgB = blank_data.mean(axis=1).tolist()

        header.append('maxB')
        header.append('avgB')
    else:
        # Otherwise, set the 'blanks' to an empty list
        blank_sample_ids = []
        blanks = False

    # If specific list of experimental sample IDs aren't provided, 
    # assume everything not marked blank is an experimental sample

    if len(exp_sample_ids) == 0:
        exp_sample_ids = set(biom_file.ids(axis='sample')) - set(blank_sample_ids)

    sample_data = biom_file.filter(exp_sample_ids, axis='sample', 
        invert=False, inplace=False).matrix_data

    maxS = sample_data.max(axis=1).todense().tolist()
    avgS = sample_data.mean(axis=1).tolist()

    stats_dict = {}

    i = 0
    if blanks:
        for otu in biom_file.ids(axis='observation'):
            stats_dict[otu] = [maxS[i][0], avgS[i][0], maxB[i][0], avgB[i][0]]
            i += 1
    else:
        for otu in biom_file.ids(axis='observation'):
            stats_dict[otu] = [maxS[i][0], avgS[i][0]]
            i += 1

    return(header, stats_dict)

def pick_min_relabund_threshold(stats_dict, stats_header, min_relabund, sample_stat='maxS'):

    i_s = stats_header.index(sample_stat)

    passed_otus = set()

    for otu in stats_dict:
        if(float(stats_dict[otu][i_s]) < float(min_relabund)):
            passed_otus.add(otu)

    return(passed_otus)

def compare_blank_abundances(stats_dict, stats_header,
                            sample_stat, blank_stat, scalar=1, negate=False):
    """Note that this method will default to returning sequences for which
    the criteria sample_stat > blank_stat * scalar are TRUE, i.e. non-contam
    sequences. To return contaminants (sequences that FAIL the inequality),
    set negate to True."""

    i_s = stats_header.index(sample_stat)
    i_b = stats_header.index(blank_stat)

    passed_otus = set()

    for otu in stats_dict:
        if((float(stats_dict[otu][i_s]) > (float(scalar) * float(stats_dict[otu][i_b]))) != negate):
            passed_otus.add(otu)

    # print passed_otus
    return(passed_otus)


def calc_per_category_decontam_stats(biom_obj, filter_otus):
    reads = biom_obj.filter(lambda val, id_, metadata: id_ in filter_otus,
                     axis='observation', invert=False, inplace=False).sum(axis = 'sample')
    otus = biom_obj.pa(inplace = False).filter(lambda val, id_, metadata: id_ in filter_otus,
                     axis='observation', invert=False, inplace=False).sum(axis = 'sample')

    return(reads.tolist(),otus.tolist())


def calc_per_library_decontam_stats(start_biom, output_dict):
    # calculate starting number of sequences and unique sequences per library

    steps = ['below_relabund_threshold','putative_contaminants','ever_good_seqs','reinstated_seqs','all_good_seqs']

    results_dict = {}

    results_dict['starting'] = calc_per_category_decontam_stats(start_biom, start_biom.ids(axis='observation'))
    results_header = ['starting']

    for step in steps:
        if step in output_dict:
            results_dict[step] = calc_per_category_decontam_stats(start_biom, output_dict[step])
            results_header.append(step)

    return(results_dict, results_header)


def filter_contaminated_libraries(unique_seq_biom, contaminant_otus, contam_threshold):
    # make relabund table
    norm_biom = unique_seq_biom.norm(inplace = False)

    # filter out sequences marked as contaminants
    norm_biom.filter(lambda val, id_, metadata: id_ in contaminant_otus,
                     axis='observation', invert=True, inplace=True)

    # filter out samples above threshold
    norm_biom.filter(lambda val, id_, metadata: sum(val) > contam_threshold,
                          axis='sample', invert=False, inplace=True)

    # filter contam sequences from original biom
    filtered_biom = unique_seq_biom.filter(lambda val, id_, metadata: id_ in contaminant_otus,
                     axis='observation', invert=True, inplace=False)

    # filter samples that lost too much relative to starting from original biom
    filtered_biom = filtered_biom.filter(lambda val, id_, metadata: id_ in norm_biom.ids(axis='sample'),
                     axis='sample', invert=False, inplace=True)

    return(filtered_biom)


def print_filtered_otu_map(input_otu_map_fp, output_otu_map_fp, filter_set):

    output_otu_map_f = open(output_otu_map_fp, 'w')
    for line in open(input_otu_map_fp, 'U'):
        seq_identifier = line.strip().split('\t')[0]
        # write OTU line if present in the filter set
        if seq_identifier in filter_set:
            output_otu_map_f.write(line)
    output_otu_map_f.close()
    return

def print_filtered_mothur_counts(mothur_counts_fp, output_counts_fp, filter_set):

    output_counts_f = open(output_counts_fp, 'w')

    t = 0

    for line in open(mothur_counts_fp, 'U'):
        seq_identifier = line.strip().split('\t')[0]

        # only write this line if the otu has more than n sequences (so
        # greater than n tab-separated fields including the otu identifier)
        # or if it's the header (first) line
        if seq_identifier in filter_set or t == 0:
            output_counts_f.write(line)
        t += 1

    output_counts_f.close()
    return


def print_per_library_stats(per_library_stats, per_library_stats_header, lib_ids, dropped_libs=[]):

    outline = 'Library\t'

    outline += '_reads\t'.join(per_library_stats_header) + '_reads\t'

    outline += '_otus\t'.join(per_library_stats_header) + '_otus'

    if len(dropped_libs) > 0:
        outline += '\tlibrary_discarded'
        discard = True
    else:
        discard = False

    outline += '\n'

    t = 0
    for lib in lib_ids:


        outline += lib

        for category in per_library_stats_header:
            outline += '\t' + str(int(per_library_stats[category][0][t]))

        for category in per_library_stats_header:
            outline += '\t' + str(int(per_library_stats[category][1][t]))

        if discard:
            if lib in dropped_libs:
                outline += '\tTrue'
            else:
                outline += '\tFalse'

        outline += '\n'

        t += 1

    return(outline)


def print_otu_disposition(input_seqs, output_dict, hierarchy=[]):
    outline = ''

    if hierarchy == []:
        hierarchy = ['below_relabund_threshold', 'putative_contaminants','reinstated_seqs','ever_good_seqs']

    # Subset hierarchy to levels also in output dict:

    hierarchy = [x for x in hierarchy if x in output_dict]

    # Check that the levels of the hierarchy are non-overlapping:

    for x in range(len(hierarchy) - 1):
        for y in range(x + 1,len(hierarchy)):
            if not output_dict[hierarchy[x]].isdisjoint(output_dict[hierarchy[y]]):
                print('warning: non-disjoint sets in the disposition hierarchy')

    seqs_left = set(input_seqs)

    for seq in input_seqs:
        for level in hierarchy:
            if seq in output_dict[level]:
                outline += '{0}\t{1}\n'.format(seq,level)
                break

    return(outline)


def print_filtered_seq_headers(seq_headers, output_headers_fp, filter_set):

    output_headers_f = open(output_headers_fp, 'w')

    for x in seq_headers:
        if x in filter_set:
            output_headers_f.write('{0}\n'.format(x))

    output_headers_f.close()
    return


def print_filtered_output(output_method, unfiltered_input, output_dir, output_dict, output_categories=None):
    output_fn = 'print_filtered_' + output_method

    if not output_categories:
        output_categories = output_dict.keys()

    if output_method == 'seq_headers':
        output_fn = print_filtered_seq_headers
    elif output_method == 'mothur_counts':
        output_fn = print_filtered_mothur_counts
    elif output_method == 'otu_map':
        output_fn = print_filtered_otu_map

    for category in output_categories:
        output_fn(unfiltered_input,
                  os.path.join(output_dir,
                          '{0}_{1}.txt'.format(category, output_method)),
                  output_dict[category])
    return


def print_results_file(seq_ids,
                       output_dict,
                       output_fp,
                       stats_header=None,
                       stats_dict=None,
                       corr_data_dict=None):

    output_f = open(output_fp, 'w')

    header = "SeqID"

    sorted_categories = sorted(output_dict.keys())

    for category in sorted_categories:
        header += '\t{0}'.format(category)

    if stats_header:
        for x in stats_header:
            header += '\t{0}'.format(x)

    if corr_data_dict:
        header += '\t{0}\t{1}'.format('spearman_r','spearman_p')

    output_f.write(header + '\n')

    for otu in seq_ids:
        outline = str(otu)

        for category in sorted_categories:
            outline += '\t{0}'.format(1 if otu in output_dict[category] else 0)

        if stats_header:
            t = 0
            for x in stats_header:
                outline += '\t{0:.3f}'.format(stats_dict[otu][t])
                t += 1

        if corr_data_dict:
            outline += '\t{0:.3f}\t{1:.3f}'.format(
                                           corr_data_dict[otu][0],
                                           corr_data_dict[otu][1])

        output_f.write(outline + '\n')

    return


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    otu_table_fp = opts.otu_table_fp
    mothur_counts_fp = opts.mothur_counts_fp
    mapping_fp = opts.mapping_fp
    valid_states = opts.valid_states
    blank_id_fp = opts.blank_id_fp
    contaminant_db_fp = opts.contaminant_db_fp
    contaminant_similarity = opts.contaminant_similarity
    max_correlation = opts.max_correlation
    correlate_header = opts.correlate_header
    input_fasta_fp = opts.input_fasta_fp
    otu_map_fp = opts.otu_map_fp
    output_dir = opts.output_dir
    min_relabund_threshold = opts.min_relabund_threshold
    prescreen_threshold = opts.prescreen_threshold
    removal_stat_blank = opts.removal_stat_blank
    removal_stat_sample = opts.removal_stat_sample
    removal_differential = opts.removal_differential
    reinstatement_stat_sample = opts.reinstatement_stat_sample
    reinstatement_stat_blank = opts.reinstatement_stat_blank
    reinstatement_differential = opts.reinstatement_differential
    reinstatement_sample_number = opts.reinstatement_sample_number
    reinstatement_method = opts.reinstatement_method
    write_output_seq_lists = opts.write_output_seq_lists
    write_filtered_output = opts.write_filtered_output
    drop_lib_threshold = opts.drop_lib_threshold
    write_per_seq_stats = opts.write_per_seq_stats
    write_per_library_stats = opts.write_per_library_stats
    write_per_seq_disposition = opts.write_per_seq_disposition

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

    # Make sure passed at least one of an OTU biom or mothur counts table file
    input_file_counter = 0

    if mothur_counts_fp:
        input_file_counter += 1
        unique_seq_biom = mothur_counts_to_biom(mothur_counts_fp)
        mothur_output = True
        print "mothur input"

    if otu_table_fp:
        input_file_counter += 1
        unique_seq_biom = load_table(otu_table_fp)
        mothur_output = False
        print "BIOM input"

    if input_file_counter != 1:
        option_parser.error("must provide ONLY ONE of an OTU table biom file or"
                            "mothur counts table")

    # Check to make sure that if blank-based contamination filtering requested,
    # all necessary options are specified:

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
    elif removal_options_counter == 0:
        blank_stats_removal = False
    elif removal_options_counter == 3:
        blank_stats_removal = True


    # If reference-based filtering requested, make sure all necessary options
    # have been specified:

    if contaminant_db_fp and not input_fasta_fp:
        option_parser.error("If specifying ref-based contaminant ID, must "
                            "also specify path to input sequence fasta")


    # If correlation-based filtering requested, make sure correlate data 
    # are specified

    if max_correlation and not correlate_header:
        option_parser.error("If specifying maximum Spearman correlation, must "
                           "also provide map column header for correlate data")


    # If sequence reinstatement is requested, make sure all necessary options
    # are specified

    reinstatement_options_counter = 0
    if reinstatement_stat_blank:
        reinstatement_options_counter += 1
    if reinstatement_stat_sample:
        reinstatement_options_counter += 1
    if reinstatement_differential:
        reinstatement_options_counter += 1

    if ((reinstatement_options_counter > 0) and 
        (reinstatement_options_counter < 3)):
        option_parser.error("Must provide all of "
                            "reinstatement_stats_blank, "
                            "reinstatement_stat_sample, and "
                            "reinstatement_differential, or none.")

    if ((reinstatement_options_counter == 3 and reinstatement_sample_number)
        and not reinstatement_method):
        option_parser.error("If providing sample number AND abundance criteria "
                            "for sequence reinstatement, must also provide "
                            "a method for combining results.")

    if reinstatement_options_counter == 3 or reinstatement_sample_number:
        reinstatement = True
    else:
        reinstatement = False

    # get blank sample IDs from mapping file or sample ID list

    if mapping_fp and valid_states:
        blank_sample_ids = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)
        blanks = True
    elif blank_id_fp is not None:
        blank_id_f = open(blank_id_fp, 'Ur')
        blank_sample_ids = set([line.strip().split()[0]
                                for line in blank_id_f
                                if not line.startswith('#')])
        blank_id_f.close()
        blanks = True
    else:
        blanks = False


    # Initialize output objets  

    output_dict = {}
    contaminant_types = []

    contamination_stats_dict = None
    contamination_stats_header = None
    corr_data_dict = None

    # Do blank-based stats calculations, if not there check to make sure no 
    # blank-dependent methods are requested:

    if blanks:
        if prescreen_threshold:
            low_contam_libraries = prescreen_libraries(unique_seq_biom,
                                                       blank_sample_ids,
                                                       removal_stat_sample, 
                                                       removal_stat_blank, 
                                                       removal_differential, 
                                                       prescreen_threshold)

            contamination_stats_header, contamination_stats_dict = \
                get_contamination_stats(unique_seq_biom,
                                        blank_sample_ids,
                                        exp_sample_ids=low_contam_libraries)
        else:
            contamination_stats_header, contamination_stats_dict = \
                get_contamination_stats(unique_seq_biom, blank_sample_ids)

    elif (blank_stats_removal or reinstatement or prescreen_threshold):
        option_parser.error("Blank-based filtering requested but no blank"
                            "samples indicated in mapping file or ID file.")
    else:
        contamination_stats_header, contamination_stats_dict = \
            get_contamination_stats(unique_seq_biom)


    seq_ids = unique_seq_biom.ids(axis='observation')


    # Do blank-based contaminant identification

    if min_relabund_threshold:
        output_dict['below_relabund_threshold'] = pick_min_relabund_threshold(
                                                  contamination_stats_dict,
                                                  contamination_stats_header,
                                                  min_relabund_threshold)


    if blank_stats_removal:
        output_dict['abund_contaminants'] = compare_blank_abundances(contamination_stats_dict, 
                                contamination_stats_header,
                                removal_stat_sample,
                                removal_stat_blank,
                                removal_differential,
                                negate=True)

        contaminant_types.append('abund_contaminants')


    # Do reference-based contaminant identification

    if contaminant_db_fp:
        output_dict['ref_contaminants'] = pick_ref_contaminants(seq_ids, contaminant_db_fp, input_fasta_fp, contaminant_similarity, output_dir)

        contaminant_types.append('ref_contaminants')


    # Do spearman correlation based contaminant identification

    if max_correlation:
        metadata_dict = parse_mapping_file_to_dict(open(mapping_fp, 'U'))[0]

        corr_data_dict = {x: float(metadata_dict[x][correlate_header]) for x in metadata_dict}

        output_dict['corr_contaminants'], corr_contaminant_dict = pick_corr_contaminants(unique_seq_biom,
                                                   corr_data_dict,
                                                   max_correlation)

        contaminant_types.append('corr_contaminants')
    else:
        corr_contaminant_dict = None


    # Putative contaminants are those that have been identified by any method

    output_dict['putative_contaminants'] = set.union(*map(set, [output_dict[x] for x in contaminant_types]))


    # If considering low abundance sequences, remove those from consideration as potential contaminants 

    if 'below_relabund_threshold' in output_dict:
        output_dict['putative_contaminants'] = output_dict['putative_contaminants'] - set(output_dict['below_relabund_threshold'])


    # Pick abundance-criterion seqs to reinstate

    if (reinstatement_stat_blank and reinstatement_stat_sample and reinstatement_differential):
        output_dict['abund_reinstated_seqs'] = reinstate_abund_seqs(output_dict['putative_contaminants'], 
                     contamination_stats_dict, 
                     contamination_stats_header,
                     reinstatement_stat_sample,
                     reinstatement_stat_blank,
                     reinstatement_differential)

        output_dict['reinstated_seqs'] = output_dict['abund_reinstated_seqs']


    # Pick incidence-criterion seqs to reinstate
    if reinstatement_sample_number:
        output_dict['incidence_reinstated_seqs'] = reinstate_incidence_seqs(
                     output_dict['putative_contaminants'],
                     unique_seq_biom,
                     blank_sample_ids,
                     reinstatement_sample_number)

        output_dict['reinstated_seqs'] = output_dict['incidence_reinstated_seqs']


    # combine incidence and abundance reinstatements
    if reinstatement_sample_number and reinstatement_stat_blank:
        if reinstatement_method == "union":
            output_dict['reinstated_seqs'] = output_dict['abund_reinstated_seqs'] | output_dict['incidence_reinstated_seqs']
        elif reinstatement_method == "intersection":
            output_dict['reinstated_seqs'] = output_dict['abund_reinstated_seqs'] & output_dict['incidence_reinstated_seqs']


    # make sets for sequence _never_ identified as contaminants:

    output_dict['ever_good_seqs'] = set(seq_ids) - output_dict['putative_contaminants']

    # If considering low abundance sequences, remove those from consideration as potential contaminants 

    if 'below_relabund_threshold' in output_dict:
        output_dict['ever_good_seqs'] = output_dict['ever_good_seqs'] - set(output_dict['below_relabund_threshold'])

    # Make set of good seqs for final filtering

    final_good_seqs = output_dict['ever_good_seqs']

    # ...and those either never ID'd as contaminants or reinstated:
    if reinstatement:
        output_dict['all_good_seqs'] = set(output_dict['ever_good_seqs'] | output_dict['reinstated_seqs'])
        final_good_seqs = output_dict['all_good_seqs']
        # ...and those who remain contaminants after reinstatement:
        output_dict['never_good_seqs'] = set(output_dict['putative_contaminants'] - output_dict['reinstated_seqs'])


    # print filtered OTU maps if given a QIIME OTU map input

    if otu_map_fp:
        print_filtered_output('otu_map', otu_map_fp, output_dir, output_dict)


    # print filtered Mothur counts tables if given a Mothur counts table input

    if mothur_output:
        print_filtered_output('mothur_counts', mothur_counts_fp, output_dir, output_dict)


    # print filtered seq header files if requested

    if write_output_seq_lists:
        print_filtered_output('seq_headers', seq_ids, output_dir, output_dict)


    # filter final biom file to just good seqs

    filtered_biom = unique_seq_biom.filter(lambda val, id_, metadata: id_ in final_good_seqs,
                     axis='observation', invert=False, inplace=False)

    # drop heavily contaminated libraries if requested

    if drop_lib_threshold:
        dropped_libs = unique_seq_biom.norm(inplace=False).filter(lambda val, id_, metadata: id_ in final_good_seqs,
                 axis='observation', invert=False, inplace=False).filter(lambda val, id_, metadata: sum(val) >= drop_lib_threshold,
                 axis='sample', invert=True, inplace=False).ids(axis='sample')
        filtered_biom.filter(lambda val, id_, metadata: id_ in dropped_libs,
                 axis='sample', invert=True, inplace=True)
    else:
        dropped_libs = []


    # print filtered biom/mothur_output if library filtering is requested

    if write_filtered_output:

        if mothur_output:
            output_counts_string = biom_to_mothur_counts(filtered_biom)
            with open(os.path.join(output_dir,'decontaminated_table.counts'), "w") as output_counts_file:
                output_counts_file.write(output_counts_string)
        else:
            output_biom_string = filtered_biom.to_json('Filtered by decontaminate.py')
            output_biom_string
            with open(os.path.join(output_dir,'decontaminated_otu_table.biom'), "w") as output_biom_file:
                output_biom_file.write(output_biom_string)



    # print per-library stats if requested

    if write_per_library_stats:
        per_library_stats, per_library_stats_header = calc_per_library_decontam_stats(unique_seq_biom, output_dict)
        library_stats_string = print_per_library_stats(per_library_stats, per_library_stats_header, unique_seq_biom.ids(axis='sample'), dropped_libs=dropped_libs)
        
        with open(os.path.join(output_dir,'decontamination_per_library_stats.txt'), "w") as output_stats_file:
            output_stats_file.write(library_stats_string)


    # print otu by disposition file if requested

    if write_per_seq_disposition:
        per_seq_disposition = print_otu_disposition(seq_ids, output_dict)

        with open(os.path.join(output_dir,'decontamination_per_otu_disposition.txt'), "w") as output_stats_file:
            output_stats_file.write(per_seq_disposition)


    # print log file / per-seq info
    if write_per_seq_stats:
        print_results_file(seq_ids,
                       output_dict,
                       os.path.join(output_dir,'contamination_summary.txt'),
                       contamination_stats_header,
                       contamination_stats_dict,
                       corr_contaminant_dict)





if __name__ == "__main__":
    main()




