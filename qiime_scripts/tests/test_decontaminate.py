#!/usr/bin/env python

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Dan Knights"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from unittest import TestCase, main

from tempfile import mkstemp, mkdtemp
import os

import numpy
from numpy import array
from numpy.testing import assert_almost_equal, assert_allclose
from StringIO import StringIO

from biom import load_table, parse_table

from decontaminate import *

def assertDeepAlmostEqual(test_case, expected, actual, *args, **kwargs):
    """
    Assert that two complex structures have almost equal contents.

    Compares lists, dicts and tuples recursively. Checks numeric values
    using test_case's :py:meth:`unittest.TestCase.assertAlmostEqual` and
    checks all other values with :py:meth:`unittest.TestCase.assertEqual`.
    Accepts additional positional and keyword arguments and pass those
    intact to assertAlmostEqual() (that's how you specify comparison
    precision).

    :param test_case: TestCase object on which we can call all of the basic
    'assert' methods.
    :type test_case: :py:class:`unittest.TestCase` object

    test from https://github.com/larsbutler/oq-engine/blob/master/tests/utils/helpers.py
    """
    is_root = not '__trace' in kwargs
    trace = kwargs.pop('__trace', 'ROOT')
    try:
        if isinstance(expected, (int, float, long, complex)):
            test_case.assertAlmostEqual(expected, actual, *args, **kwargs)
        elif isinstance(expected, (list, tuple, numpy.ndarray)):
            test_case.assertEqual(len(expected), len(actual))
            for index in xrange(len(expected)):
                v1, v2 = expected[index], actual[index]
                assertDeepAlmostEqual(test_case, v1, v2,
                                      __trace=repr(index), *args, **kwargs)
        elif isinstance(expected, dict):
            test_case.assertEqual(set(expected), set(actual))
            for key in expected:
                assertDeepAlmostEqual(test_case, expected[key], actual[key],
                                      __trace=repr(key), *args, **kwargs)
        else:
            test_case.assertEqual(expected, actual)
    except AssertionError as exc:
        exc.__dict__.setdefault('traces', []).append(trace)
        if is_root:
            trace = ' -> '.join(reversed(exc.traces))
            exc = AssertionError("%s\nTRACE: %s" % (exc.message, trace))
        raise exc


class DecontaminationTests(TestCase):

    def setUp(self):
        """Init variables for the tests """

        self.test_biom = parse_table(test_biom_file)

    def test_pick_min_relabund_threshold(self):

        contamination_header = exp_contamination_header

        contamination_stats_dict = exp_contamination_stats_dict_prop

        exp_below_threshold_seqs = set(['contam1','contam4'])

        obs_below_threshold_seqs = pick_min_relabund_threshold(contamination_stats_dict, exp_contamination_header, 0.1)

        self.assertEqual(exp_below_threshold_seqs,obs_below_threshold_seqs)

    def test_pick_corr_contaminants(self):

        corr_data_dict = {'Blank1': 0, 'Blank2': 1, 'Sample1': 2, 'Sample2': 3}

        test_biom = self.test_biom

        exp_corr_contaminants = set(['contam1','contam2','contam4'])

        exp_corr_dict = {u'otu1': (0.8, 0.2),
                         u'otu2': (0.6, 0.4),
                         u'otu3': (-0.4, 0.6),
                         u'otu4': (0.2581988897471611, 0.7418011102528389),
                         u'contam1': (-0.8, 0.2),
                         u'contam2': (-1.0, 0.0),
                         u'contam3': (0.4, 0.6),
                         u'contam4': (-0.89442719099991586, 0.10557280900008413)}

        obs_corr_contaminants, obs_corr_dict = pick_corr_contaminants(test_biom,
                                                           corr_data_dict,
                                                           -0.6)

        self.assertEqual(obs_corr_contaminants, exp_corr_contaminants)
        assertDeepAlmostEqual(self, obs_corr_dict, exp_corr_dict)


    def test_reinstate_abund_seqs(self):

        putative_contaminants = set(['otu1','otu2','otu3','contam1','contam2','contam3','contam4'])

        contamination_stats_dict = exp_contamination_stats_dict

        contamination_stats_header = exp_contamination_header


        exp_abund_reinstated_seqs_max = set(['otu1','otu2','otu3'])

        obs_abund_reinstated_seqs_max = reinstate_abund_seqs(putative_contaminants, 
                         contamination_stats_dict, 
                         contamination_stats_header,
                         'maxS',
                         'maxB',
                         1)

        self.assertEqual(exp_abund_reinstated_seqs_max, obs_abund_reinstated_seqs_max)


        exp_abund_reinstated_seqs_avg = set(['otu1','otu2','contam3'])

        obs_abund_reinstated_seqs_avg = reinstate_abund_seqs(putative_contaminants, 
                         contamination_stats_dict, 
                         contamination_stats_header,
                         'avgS',
                         'avgB',
                         1)

        self.assertEqual(exp_abund_reinstated_seqs_avg, obs_abund_reinstated_seqs_avg)

    def test_reinstate_incidence_seqs(self):

        test_biom = self.test_biom

        blank_sample_ids = ['Blank1', 'Blank2']

        putative_contaminants = set(['otu1','otu2','otu3','contam1','contam2','contam3','contam4'])

        exp_incidence_reinstated_seqs = set(['otu1','otu2','contam3'])

        obs_incidence_reinstated_seqs = reinstate_incidence_seqs(
                         putative_contaminants,
                         test_biom,
                         blank_sample_ids,
                         2)

        self.assertEqual(exp_incidence_reinstated_seqs, obs_incidence_reinstated_seqs)


    def test_mothur_counts_to_biom(self):
        """Convert a mothur counts table to a biom object"""

        test_counts_f = test_mothur_counts_table

        exp_counts_biom = self.test_biom

        obs_counts_biom = mothur_counts_to_biom(StringIO(test_counts_f))

        self.assertEqual(exp_counts_biom, obs_counts_biom)


    def test_biom_to_mothur_counts(self):
        """Convert a biom OTU table object to a mothur counts table"""
        test_biom = self.test_biom

        exp_counts_table = test_mothur_counts_table

        obs_counts_table = biom_to_mothur_counts(test_biom)

        self.assertEqual(exp_counts_table,obs_counts_table)


    def test_filter_contaminated_libraries(self):
        """Test filtering process for heavily contaminated libraries"""

        test_biom = self.test_biom
        exp_filtered_biom = parse_table(test_biom_file_filt)

        contaminant_otus = set(['otu3','contam1','contam2','contam4'])

        obs_filtered_biom = filter_contaminated_libraries(test_biom, contaminant_otus, 0.8)

        self.assertEqual(exp_filtered_biom, obs_filtered_biom)


    def test_print_otu_disposition(self):

        output_dict = {'below_relabund_threshold': set(['otu4','contam4']),
                       'putative_contaminants': set(['contam1','contam2']),
                       'ever_good_seqs': set(['otu1','otu2','otu3','contam3'])}

        input_seqs = ['otu1','otu2','otu3','otu4','contam1','contam2','contam3','contam4']

        exp_output_string = 'otu1\tever_good_seqs\notu2\tever_good_seqs\notu3\tever_good_seqs\notu4\tbelow_relabund_threshold\ncontam1\tputative_contaminants\ncontam2\tputative_contaminants\ncontam3\tever_good_seqs\ncontam4\tbelow_relabund_threshold\n'

        obs_output_string = print_otu_disposition(input_seqs, output_dict)

        self.assertEqual(exp_output_string,obs_output_string)


    def test_calc_per_library_decontam_stats(self):

        test_biom = self.test_biom

        output_dict = {'below_relabund_threshold': set(['otu4','contam4']),
                       'putative_contaminants': set(['contam1','contam2']),
                       'ever_good_seqs': set(['otu1','otu2','otu3','contam3'])}

        exp_results_dict = {'starting': ([24.0, 8.0, 26.0, 5.0], [7.0, 5.0, 7.0, 3.0]),
                            'below_relabund_threshold': ([3.0, 1.0, 3.0, 0.0], [1.0, 1.0, 1.0, 0.0]),
                            'putative_contaminants': ([8.0, 3.0, 4.0, 0.0], [2.0, 2.0, 2.0, 0.0]),
                            'ever_good_seqs': ([13.0, 4.0, 19.0, 5.0], [4.0, 2.0, 4.0, 3.0])}

        exp_results_header = ['starting','below_relabund_threshold', 'putative_contaminants', 'ever_good_seqs']

        obs_results_dict, obs_results_header = calc_per_library_decontam_stats(test_biom, output_dict)

        assertDeepAlmostEqual(self, exp_results_dict, obs_results_dict)
        self.assertEqual(exp_results_header, obs_results_header)


    def test_prescreen_libraries(self):
        """Test pre-screening libraries for contamination"""

        unique_seq_biom = self.test_biom

        exp_exclude_libs = set(['Sample2'])

        obs_exclude_libs = prescreen_libraries(unique_seq_biom,
                                               blank_sample_ids = ['Blank1','Blank2'],
                                               removal_stat_sample = 'maxS',
                                               removal_stat_blank = 'maxB',
                                               removal_differential = 1,
                                               prescreen_threshold = 0.65)

        self.assertEqual(exp_exclude_libs, set(obs_exclude_libs))


        exp_exclude_libs = set(['Sample1', 'Sample2'])

        obs_exclude_libs = prescreen_libraries(unique_seq_biom,
                                               blank_sample_ids = ['Blank1','Blank2'],
                                               removal_stat_sample = 'maxS',
                                               removal_stat_blank = 'maxB',
                                               removal_differential = 1,
                                               prescreen_threshold = 0.5)

        self.assertEqual(exp_exclude_libs, set(obs_exclude_libs))


    def test_pick_ref_contaminants(self):
        """Test the reference-based contaminant search"""

        # make temp input and output files 

        test_output_dir = mkdtemp()


        _, tmp_ref_fasta_filepath = mkstemp(prefix="decontaminate_ref",
                                                    suffix=".fasta")

        tmp_ref_fasta = open(tmp_ref_fasta_filepath, "w")
        tmp_ref_fasta.write(blanks_rep_set_fasta)
        tmp_ref_fasta.close()

        _, tmp_test_fasta_filepath = mkstemp(prefix="decontaminate_test",
                                                    suffix=".fasta")

        tmp_test_fasta = open(tmp_test_fasta_filepath, "w")
        tmp_test_fasta.write(unique_seqs_rep_set_fasta)
        tmp_test_fasta.close()


        test_queries = set(['denovo47', 'denovo46', 'denovo27', 'denovo26', 'denovo25', 'denovo24', 'denovo23', 'denovo22', 'denovo21', 'denovo20', 'denovo49', 'denovo48', 'denovo29', 'denovo28', 'denovo0', 'denovo1', 'denovo2', 'denovo3', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo8', 'denovo9', 'denovo18', 'denovo19', 'denovo12', 'denovo13', 'denovo10', 'denovo11', 'denovo16', 'denovo17', 'denovo14', 'denovo15', 'denovo34', 'denovo35', 'denovo36', 'denovo37', 'denovo30', 'denovo31', 'denovo32', 'denovo33', 'denovo38', 'denovo39', 'denovo50', 'denovo51', 'denovo41', 'denovo40', 'denovo43', 'denovo42', 'denovo52', 'denovo45', 'denovo44'])

        # Test at 97% ID level

        exp_hits_97 = set(['denovo20', 'denovo41', 'denovo26', 'denovo25', 'denovo24', 'denovo23', 'denovo44', 'denovo46', 'denovo48', 'denovo28', 'denovo0', 'denovo1', 'denovo2', 'denovo3', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo8', 'denovo9', 'denovo18', 'denovo19', 'denovo12', 'denovo13', 'denovo10', 'denovo16', 'denovo17', 'denovo14', 'denovo15', 'denovo34', 'denovo35', 'denovo36', 'denovo30', 'denovo31', 'denovo33', 'denovo38', 'denovo39', 'denovo51', 'denovo40', 'denovo42', 'denovo52', 'denovo45'])

        obs_hits_97 = pick_ref_contaminants(test_queries, tmp_ref_fasta_filepath, tmp_test_fasta_filepath, 0.97, test_output_dir)

        self.assertEqual(exp_hits_97, obs_hits_97)

        # Test at 99% ID level

        exp_hits_99 = set(['denovo20', 'denovo41', 'denovo40', 'denovo25', 'denovo24', 'denovo23', 'denovo46', 'denovo48', 'denovo0', 'denovo2', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo12', 'denovo13', 'denovo10', 'denovo14', 'denovo15', 'denovo34', 'denovo36', 'denovo33', 'denovo52', 'denovo39', 'denovo51', 'denovo42'])

        obs_hits_99 = pick_ref_contaminants(test_queries, tmp_ref_fasta_filepath, tmp_test_fasta_filepath, 0.99, test_output_dir)

        self.assertEqual(exp_hits_99, obs_hits_99)

        # Test at 100% ID level

        exp_hits_100 = set(['denovo14', 'denovo25', 'denovo42', 'denovo20', 'denovo48'])

        obs_hits_100 = pick_ref_contaminants(test_queries, tmp_ref_fasta_filepath, tmp_test_fasta_filepath, 1.00, test_output_dir)

        self.assertEqual(exp_hits_100, obs_hits_100)

        # Remove temporary files
        os.removedirs(test_output_dir)
        os.remove(tmp_ref_fasta_filepath)
        os.remove(tmp_test_fasta_filepath)


    def test_get_contamination_stats(self):
        """add_contamination_stats_to_biom: 
        """

        blank_sample_ids = ['Blank1', 'Blank2']

        test_biom = self.test_biom


        # test when passing already-proportional table

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, blank_sample_ids, proportional=True)

        # Header is as expected
        self.assertItemsEqual(exp_contamination_header,
                              obs_contamination_header)

        # Contamination stats dict is as expected
        self.assertDictEqual(exp_contamination_stats_dict,
                              obs_contamination_stats_dict)
        

        # test when passing non-proportional table

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, blank_sample_ids, proportional=False)

        # Header is as expected
        self.assertItemsEqual(exp_contamination_header,
                              obs_contamination_header)

        # Contamination stats dict is as expected
        assertDeepAlmostEqual(self,exp_contamination_stats_dict_prop,
                              obs_contamination_stats_dict)


        # test when passing no blank information

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, proportional=True)

        # Header is as expected
        self.assertItemsEqual(['maxS','avgS'],
                              obs_contamination_header)

        # Contamination stats dict is as expected
        assertDeepAlmostEqual(self, exp_contam_stats_dict_noblanks,
                              obs_contamination_stats_dict)


        # test when passing specified samples information

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, blank_sample_ids, exp_sample_ids=['Sample2'], proportional=True)

        # Header is as expected
        self.assertItemsEqual(exp_contamination_header,
                              obs_contamination_header)

        # Contamination stats dict is as expected
        assertDeepAlmostEqual(self, exp_contam_stats_dict_sample,
                              obs_contamination_stats_dict)


    def test_compare_blank_abundances(self):
        """examine contamination stats dict for sequences that fail a specificied
        abundance threshold for samples vs blank abundances.
        """

        contamination_stats_dict = exp_contamination_stats_dict

        contamination_header = exp_contamination_header

        # test for maxS > maxB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'maxB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu3','otu4'],passed_seqs)

        # test for maxS > maxB w/ negate
        rejected_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'maxB',
                                 scalar = 1,
                                 negate = True)
        
        self.assertItemsEqual(['contam1','contam2','contam3','contam4'],rejected_seqs)

        # test for maxS > avgB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'avgB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu3','otu4','contam2','contam3'],passed_seqs)

        # test for maxS > avgB w/ negate
        rejected_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'avgB',
                                 scalar = 1,
                                 negate = True)
        
        self.assertItemsEqual(['contam1','contam4'],rejected_seqs)

        # test for avgS > avgB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'avgS',
                                 blank_stat = 'avgB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu4','contam3'],passed_seqs)


        # test for avgS > maxB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 contamination_header,
                                 sample_stat = 'avgS',
                                 blank_stat = 'maxB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu4'],passed_seqs)


test_biom_file = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","matrix_type": "sparse","generated_by": "BIOM-Format 2.0.1-dev","date": "2014-09-14T17:26:13.141629","type": "OTU table","matrix_element_type": "float","shape": [8, 4],"data": [[0,0,1.0],[0,2,4.0],[0,3,2.0],[1,0,4.0],[1,1,1.0],[1,2,6.0],[1,3,1.0],[2,0,4.0],[2,1,3.0],[2,2,6.0],[3,2,3.0],[4,0,4.0],[4,1,2.0],[4,2,1.0],[5,0,4.0],[5,1,1.0],[5,2,3.0],[6,0,4.0],[6,2,3.0],[6,3,2.0],[7,0,3.0],[7,1,1.0]],"rows": [{"id": "otu1", "metadata": null},{"id": "otu2", "metadata": null},{"id": "otu3", "metadata": null},{"id": "otu4", "metadata": null},{"id": "contam1", "metadata": null},{"id": "contam2", "metadata": null},{"id": "contam3", "metadata": null},{"id": "contam4", "metadata": null}],"columns": [{"id": "Blank1", "metadata": null},{"id": "Blank2", "metadata": null},{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null}]}"""
test_biom_file_filt = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","matrix_type": "sparse","generated_by": "BIOM-Format 2.0.1-dev","date": "2015-07-13T11:35:24.327048","type": "OTU table","matrix_element_type": "float","shape": [4, 1],"data": [[0,0,2.0],[1,0,1.0],[3,0,2.0]],"rows": [{"id": "otu1", "metadata": null},{"id": "otu2", "metadata": null},{"id": "otu4", "metadata": null},{"id": "contam3", "metadata": null}],"columns": [{"id": "Sample2", "metadata": null}]}"""


test_mothur_counts_table = "Representative_Sequence\ttotal\tBlank1\tBlank2\tSample1\tSample2\notu1\t7\t1\t0\t4\t2\notu2\t12\t4\t1\t6\t1\notu3\t13\t4\t3\t6\t0\notu4\t3\t0\t0\t3\t0\ncontam1\t7\t4\t2\t1\t0\ncontam2\t8\t4\t1\t3\t0\ncontam3\t9\t4\t0\t3\t2\ncontam4\t4\t3\t1\t0\t0\n"
exp_contamination_header = ['maxS', 'avgS','maxB','avgB']
exp_contamination_stats_dict = {'otu1': [4,3,1,0.5], 'otu2': [6,3.5,4,2.5], 'otu3': [6,3,4,3.5], 'otu4': [3,1.5,0,0], 'contam1': [1,0.5,4,3], 'contam2': [3,1.5,4,2.5], 'contam3': [3,2.5,4,2], 'contam4': [0,0,3,2]}
exp_contam_stats_dict_noblanks = {'otu1': [4.00,1.75], 'otu2': [6.00,3.00], 'otu3': [6.00,3.25], 'otu4': [3.00,0.75], 'contam1': [4.00,1.75], 'contam2': [4.00,2.00], 'contam3': [4.00,2.25], 'contam4': [3.00,1.00]}
exp_contam_stats_dict_sample = {'otu1': [2.0,2.0, 1.0, 0.5],'otu2': [1.0,1.0, 4.0, 2.5],'otu3': [0.0,0.0, 4.0, 3.5],'otu4': [0.0,0.0, 0.0, 0.0],'contam1': [0.0,0.0, 4.0, 3.0],'contam2': [0.0,0.0, 4.0, 2.5],'contam3': [2.0,2.0, 4.0, 2.0],'contam4': [0.0,0.0, 3.0, 2.0]}
exp_contamination_stats_dict_prop = {'otu1': [0.4,0.276923076923077,0.0416666666666667,0.0208333333333333], 'otu2': [0.230769230769231,0.215384615384615,0.166666666666667,0.145833333333333], 'otu3': [0.230769230769231,0.115384615384615,0.375,0.270833333333333], 'otu4': [0.115384615384615,0.0576923076923077,0,0], 'contam1': [0.0384615384615385,0.0192307692307692,0.25,0.208333333333333], 'contam2': [0.115384615384615,0.0576923076923077,0.166666666666667,0.145833333333333], 'contam3': [0.4,0.257692307692308,0.166666666666667,0.0833333333333333], 'contam4': [0,0,0.125,0.125]} 

test_sample_map = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tBlank\tAbund\tDescription\nBlank1\tAACCGCAT\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t1\t0\tBlank1\nBlank2\tAACCGCAC\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t1\t1\tBlank2\nSample1\tAACCGCAA\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t0\t2\tSample1\nSample2\tAACCGCAG\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t0\t3\tSample2\n"""

blanks_rep_set_fasta = """>denovo0 SHNZ_91
ACGAACGCTGGCTTCGTGCCCAACACATGCAAGTCGTTCGGAAAGGCCCTGCCCCCGTAAAATGCTCGAGTGGCGAACCCCTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo1 SHNO_84
GGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCTTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo2 SHNO_93
GGACGAACGCTCCCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCAAACGGGGTACACGAGCTTCGAACGGGTGAGTAACACGTGAGCAATCTGAAATTCACTCTCCCATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo3 SHNO_89
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo4 SHNZ_61
ACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
"""

unique_seqs_rep_set_fasta = """>denovo0 SHOA_53
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGTTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo1 SHOA_52
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGGCTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo10 SHOG_21
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAGGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo11 SHOG_22
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAATTTTAACGGAATACCTTCGGGTAGGAAGATAAAAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGCAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCCATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCG
>denovo12 SHOG_48
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGGTGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo13 SHNO_92
GGACGAACGCTGGCGGCGTGCTTTTCACATGCAAGTCGAGCGGTAAGGCCTTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo14 SHNO_93
GGACGAACGCTCCCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCAAACGGGGTACACGAGCTTCGAACGGGTGAGTAACACGTGAGCAATCTGAAATTCACTCTCCCATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo15 SHOG_40
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCATGTCTGATGGAGCAACGCCGCGTG
>denovo16 SHOG_41
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAGAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAAGGGCCTACCAAGGCGATGATGCATAGCCGAGTTGAGAGACTGACCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo17 SHOG_47
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGGCTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCGACGCCGCG
>denovo18 SHOG_44
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAATTTTGACGGAATACTTCGGTAGGAAGTCAGAAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo19 SHOA_56
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAATCTTAACGGAATACTTCGGTAGGAAGTCAAGAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo2 SHOA_50
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGGCTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo20 SHNZ_91
ACGAACGCTGGCTTCGTGCCCAACACATGCAAGTCGTTCGGAAAGGCCCTGCCCCCGTAAAATGCTCGAGTGGCGAACCCCTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo21 SHNX_73
GGATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGGGGTGCTTGCACCCAGTGGCGAACGGGTGAGTAACACGTATCTAATCTACCCATTAGCGGGGGATAACAGTTGGAAACGACTGCTAATACCGCATACGACATTTTCTGGCATCAGAGAATGTTAAAAGGTTCGTTTGGATCACTAATGGATGAAGATGCGGCGCATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATACGTAGCCGAACTGAGAGGTTGATCGGCCACATCGGGACTGAGACACGGCCCGAACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGACTGAAGACGGTCTTCGGATTGTAAAGGTC
>denovo22 SHNX_72
GGATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGGGGTGCTTGCACCCAGTGGCGAACGGGTGAGTAACACGTATCTAATCTACCCATTAGCGGGGGATAACAGTTGGAAACGACTGCTAATACCGCATACGACATTTTCTGGCATCAGAGAATGTTAAAAGGTTCGTTTGGATCACTAATGGATGAAGATGCGGCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATACGTAGCCGAACTGAGGGGTTGATCGGCCACATCGGGACTGAGACACGGCCCGAACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGACTGAAGACGGTCTTCGGATTGTAAAGGTC
>denovo23 SHOG_37
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGGGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo24 SHOA_60
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTAGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo25 SHNZ_61
ACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo26 SHOG_19
GGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCGATGATGCATAGCCGAGTTGAGAGACTGACCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGC
>denovo27 SHOG_18
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGGGAAGGACATGAATTTTTCGGAAGGATTGTTTAGACCGAGCGGCGGATGGGTGAGTAACACGTAGGGAACCTGCCAAACAGACGGGGATACCACTTGGAAACAAGTGCTAATACCGGATAGAGCACTTTATCGCATGATAGAGTGAGGAAAGGGCGGCGAAAGCTGTCGCTGATTGATGGACCTGCGGCGTATTAGCTAGTTGGGGAGGTAAAGGCTCACCAAGGCGATGATACGTAGCCGACCTGAGAGGGTAAACGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATTTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGT
>denovo28 SHOG_27
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGCAAGGTAACGGCCTACCAAGGCGATGATGCATAGCCGAGTTGAGAGACTGACCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo29 SHNX_80
GGATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGGGGTGCTTGCACCCAGTGGCGAACGGGTGAGTAACACGTATCTAATCTACCCATTAGCGGGGGATAACAGTTGGAAACGACTGCTAATACCGCATACGACATTTTCTGGCATCAGAGAATGTTAAAAGGTTCGTTTGGATCACTAATGGATGAAGATGCGGCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATACGTAGCCGAACTGAGAGGTTGATCGGCCACATCGGGACTGAGACACGGCCCGAACTTCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGACTGAAGACGGTCTTCGGATTGTAAAGGTC
>denovo3 SHOA_57
ACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGAGTG
>denovo30 SHOG_38
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGGCGCAAGTCTGATGGAGCAACGCCGCG
>denovo31 SHOG_11
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGTTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCGATGATGCATAGCCGAGTTGAGAGACTGACCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo32 SHOG_10
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAATTTTAACGGAATACCTTCGGGTAGGAAGATAAAAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGCAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo33 SHOG_13
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGCGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo34 SHOG_34
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCTGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGA
>denovo35 SHOG_15
GGACGGACGCTGGCGGCGTGCGTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCGTGGGATACTTATTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGC
>denovo36 SHOG_14
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAATAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo37 SHOG_17
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAATTTTAACGGAATACCTTCGGGTAGGAAGATAAAAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGCAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo38 SHOG_16
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATAATTATTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGC
>denovo39 SHNO_87
ACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGAGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo4 SHOG_29
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCTATTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo40 SHNO_86
ACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCTTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCATGCCAGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACCTC
>denovo41 SHNO_85
GGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCCTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo42 SHNO_84
GGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGTAAGGCCTTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo43 SHNX_68
GGATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGGGGTGCTTGCACCCAGTGGCGAACGGGTGAGTAACACGTATCTAATCTACCCATTAGCGGGGGATAACAGTTGGAAACGACTGCTAATACCGCATACGACATTTTCTGGCATCAGAGAATGTTAAAAGGTTCGTTTGGATCACTAATGGATGAAGATGCGGCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATACGTAGCCGAACTGAGAGGTTGATCGGCCACATCGGGACTGAGACACGGCCCGAACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGACTGAAGACGGTCTTCGGATTGTAAAGGTC
>denovo44 SHOG_9
CTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGTGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGC
>denovo45 SHOG_5
GAACGCTGGCGGCGTGCCTAACACATGCAAGTCGGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGA
>denovo46 SHOG_7
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo47 SHOG_6
ATGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGGATCCAGGCAGCTTGCTGTCTGGTGAGAGTGGCGAACGGGTGAGTAATGCGTGACCAACCTGCCCCATACTCCGGAATAGCTCCTGGAAACGGGTGGTAATGCCGGATGCTCCGCATCATCGCATGATGGTGCGGGAAAGGGTTTACCGGCATGGGATGGGGTCGCGTCCTATCAGCTTGTTGGCGGGGTGATGGCCTGCCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo48 SHOG_1
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTACCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGAGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo49 SHNX_65
ATGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGGGGTGCTTGCACCCAGTGGCGAACGGGTGAGTAACACGTATCTAATCTACCCATTAGCGGGGGATAACAGTTGGAAACGACTGCTAATACCGCATACGACATTTTCTGGCATCAGAGAATGTTAAAAGGTTCGTTTGGATCACTAATGGATGAAGATGCGGCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATACGTAGCCGAACTGAGAGGTTGATCGGCCACATCGGGACTGAGACACGGCCCGAACTCCTACGGGAGGCAGCAGTAGGGAATTTTTCACAATGGGCGAAAGCCTGATGGAGCAATGCCGCGTGACTGAAGACGGTCTTCGGATTGTAAAGGTCTG
>denovo5 SHOA_55
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCACATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo50 SHOG_2
GGATGAAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGGATCCAGGCAGCTTGCTGTCTGGTGAGAGTGGCGAACGGGTGAGTAATGCGTGACCAACCTGCCCCATACTCCGGAATAGCTCCTGGAAACGGGTGGTAATGCCGGGTGTTCCGCATCATCGCATGATGGTGTGGGAAAGGGTTTACCGGTATGGGATGGGGTCGCGTCCTATCAGCTTGTTGGTGGGGTGATGGCCTGCCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGCGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGCGGGATGACGGCCTTCGGGTT
>denovo51 SHNZ_90
ACGAACGCTGGCAACGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTA
>denovo52 SHOG_97
GGACGAACGCTGGCGGCGTGCTTACCACATGCAAGTCGAGCGGTAAGGCCTTTCGGGGTACACGAGCGGCGAACGGGTGAGTAACACGTGAGCAATCTGCCCTTCACTCTGGGATAGCCACCGGAAACGGTGATTAATACCGGATATGACCACGCCTGGCATCTGGTGTGGTGGAAAGTTCTGGCGGTGAAGGATGGGCTCGCGGCCTATCAGCTTGTTGGTGAGGTAATGGCTCACCAAGGCTTCGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCAACGCCGCGTGAGGGATGACGGCCTTCGGGTTGTAAACC
>denovo6 SHOA_54
ACGAGCGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo7 SHOA_59
ACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAATCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCGTG
>denovo8 SHOA_58
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAACAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGTGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
>denovo9 SHOG_20
GGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAGCAAAGTTAAAGGAATACTTCGGTAGGAATTTAATAGCGCGAGCGGCGGATGGGTGAGTAACACGTGGGCAACCTGCCCTTTAGCTTGGGATACCACTTGGAAATAGGTGCTAATACCAAATAAGAAGTAAGAGCGCATGCTCAAGCTATGAAAGGCGGCTTTCGAGCTGTCACTAAAGGATGGGCCCGCGGTGCATTAGCTAGTTGGTAAGGTAACGGCTTACCAAGGCAATGATGCATAGCCGAGTTGAGAGACTGATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGCAAGTCTGATGGAGCAACGCCGCG
"""

# run tests if called from command line
if __name__ == "__main__":
    main()