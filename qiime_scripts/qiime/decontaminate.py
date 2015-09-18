#!/usr/bin/env python
#file decontaminate.py: helper functions for removing contaminants

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"

from biom import load_table, parse_table
from bfillings.uclust import get_clusters_from_fasta_filepath
from bfillings.usearch import usearch_qf
from scipy.stats import spearmanr
import os.path
import numpy as np

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
        save_uc_files=False,
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

def mothur_counts_to_biom(mothur_f):

    mothur_biom = parse_table(mothur_f)
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
            get_contamination_stats(unique_seq_biom, blank_sample_ids=blank_sample_ids)

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


def get_contamination_stats(biom_file, qS=50, qB=50, interpolation='midpoint', blank_sample_ids=None, exp_sample_ids=[], proportional=False):
    if not proportional:
        biom_file = biom_file.norm(inplace=False)

    header = ['maxS','avgS','q%sS' % qS]

    # Calculate blank stats if blank sample names are provided
    if blank_sample_ids:
        blanks = True

        blank_sample_ids = set(blank_sample_ids) & set(biom_file.ids(axis='sample'))

        blank_data = biom_file.filter(blank_sample_ids, axis='sample', 
                                      invert=False, inplace=False).matrix_data
        maxB = [x[0] for x in blank_data.max(axis=1).todense().tolist()]
        avgB = [x[0] for x in blank_data.mean(axis=1).tolist()]
        quantB = np.percentile(blank_data.todense(),qB,axis=1, interpolation=interpolation).tolist()

        header.append('maxB')
        header.append('avgB')
        header.append('q%sB' % qB)

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

    maxS = [x[0] for x in sample_data.max(axis=1).todense().tolist()]
    avgS = [x[0] for x in sample_data.mean(axis=1).tolist()]
    quantS = np.percentile(sample_data.todense(),qS,axis=1, interpolation=interpolation).tolist()
    
    stats_dict = {}

    i = 0
    if blanks:
        for otu in biom_file.ids(axis='observation'):
            stats_dict[otu] = [maxS[i], avgS[i], quantS[i], maxB[i], avgB[i], quantB[i]]
            i += 1
    else:
        for otu in biom_file.ids(axis='observation'):
            stats_dict[otu] = [maxS[i], avgS[i], quantS[i]]
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
                outline += '\t{0:.10f}'.format(stats_dict[otu][t])
                t += 1

        if corr_data_dict:
            outline += '\t{0:.3f}\t{1:.3f}'.format(
                                           corr_data_dict[otu][0],
                                           corr_data_dict[otu][1])

        output_f.write(outline + '\n')

    return


