#!/usr/bin/env python
#file decontaminate.py: helper functions for removing contaminants

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"

from biom import load_table
from bfillings.uclust import get_clusters_from_fasta_filepath
from bfillings.usearch import usearch_qf
from scipy.stats import spearmanr
import os.path

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

def get_contamination_stats(biom_file, blank_sample_ids=None, exp_sample_ids=None, proportional=False):
    if not proportional:
        biom_file = biom_file.norm()

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

    if not exp_sample_ids:
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
    
    i_s = stats_header.index(sample_stat)
    i_b = stats_header.index(blank_stat)

    passed_otus = set()

    for otu in stats_dict:
        if(float(stats_dict[otu][i_s]) > (float(scalar) * float(stats_dict[otu][i_b]))):
            passed_otus.add(otu)

    # print passed_otus
    return(passed_otus)

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
        # print header line
        if t == 0:
            t += 1
            continue

        seq_identifier = line.strip().split('\t')[0]

        # only write this line if the otu has more than n sequences (so
        # greater than n tab-separated fields including the otu identifier)
        if seq_identifier in filter_set:
            output_counts_f.write(line)
        t += 1

    output_counts_f.close()
    return


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
                outline += '\t{0}'.format(stats_dict[otu][t])
                t += 1

        if corr_data_dict:
            outline += '\t{0}\t{1}'.format(
                                           corr_data_dict[otu][0],
                                           corr_data_dict[otu][1])

        output_f.write(outline + '\n')

    return


