#!/usr/bin/env python
#file decontaminate.py: helper functions for removing contaminants

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"

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


