#!/usr/bin/env python

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Dan Knights"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from unittest import TestCase, main

import numpy
from numpy import array
from numpy.testing import assert_almost_equal, assert_allclose


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
        test_biom_fp = '/Users/jonsanders/Development/git_sw/qiime/qiime_test_data/decontaminate/test_otu_table.biom'

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
        test_counts_fp = '/Users/jonsanders/Development/git_sw/decontaminate/qiime_scripts/qiime_test_data/decontaminate/test.count_table'

        test_counts_biom_fp = '/Users/jonsanders/Development/git_sw/decontaminate/qiime_scripts/qiime_test_data/decontaminate/test.count_table.biom'

        test_counts_biom = load_table(test_counts_biom_fp)

        obs_counts_biom = mothur_counts_to_biom(test_counts_fp)

        self.assertEqual(test_counts_biom, obs_counts_biom)


    def test_prescreen_libraries(self):
        """Test pre-screening libraries for contamination"""

        unique_seq_biom = self.test_biom

        exp_exclude_libs = set(['Blank1', 'Blank2', 'Sample1'])

        obs_exclude_libs = prescreen_libraries(unique_seq_biom,
                                               blank_sample_ids = ['Blank1','Blank2'],
                                               removal_stat_blank = 'maxB',
                                               removal_stat_sample = 'maxS',
                                               removal_differential = 1,
                                               prescreen_threshold = 0.65)

        self.assertEqual(exp_exclude_libs, set(obs_exclude_libs))


        exp_exclude_libs = set(['Blank1', 'Blank2'])

        obs_exclude_libs = prescreen_libraries(unique_seq_biom,
                                               blank_sample_ids = ['Blank1','Blank2'],
                                               removal_stat_blank = 'maxB',
                                               removal_stat_sample = 'maxS',
                                               removal_differential = 1,
                                               prescreen_threshold = 0.5)

        self.assertEqual(exp_exclude_libs, set(obs_exclude_libs))

    def test_pick_ref_contaminants(self):
        """Test the reference-based contaminant search"""

        test_queries = set(['denovo47', 'denovo46', 'denovo27', 'denovo26', 'denovo25', 'denovo24', 'denovo23', 'denovo22', 'denovo21', 'denovo20', 'denovo49', 'denovo48', 'denovo29', 'denovo28', 'denovo0', 'denovo1', 'denovo2', 'denovo3', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo8', 'denovo9', 'denovo18', 'denovo19', 'denovo12', 'denovo13', 'denovo10', 'denovo11', 'denovo16', 'denovo17', 'denovo14', 'denovo15', 'denovo34', 'denovo35', 'denovo36', 'denovo37', 'denovo30', 'denovo31', 'denovo32', 'denovo33', 'denovo38', 'denovo39', 'denovo50', 'denovo51', 'denovo41', 'denovo40', 'denovo43', 'denovo42', 'denovo52', 'denovo45', 'denovo44'])

        test_ref_db_fp = '/Users/jonsanders/Development/git_sw/decontaminate/qiime_scripts/qiime_test_data/decontaminate/blanks_rep_set.filtered.fna'

        test_input_fasta_fp = '/Users/jonsanders/Development/git_sw/decontaminate/qiime_scripts/qiime_test_data/decontaminate/unique_seqs_rep_set.fna'

        test_output_dir = '/Users/jonsanders/Development/git_sw/decontaminate/qiime_scripts/qiime_test_data/decontaminate/testout'

        # Test at 97% ID level

        exp_hits_97 = set(['denovo20', 'denovo41', 'denovo26', 'denovo25', 'denovo24', 'denovo23', 'denovo44', 'denovo46', 'denovo48', 'denovo28', 'denovo0', 'denovo1', 'denovo2', 'denovo3', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo8', 'denovo9', 'denovo18', 'denovo19', 'denovo12', 'denovo13', 'denovo10', 'denovo16', 'denovo17', 'denovo14', 'denovo15', 'denovo34', 'denovo35', 'denovo36', 'denovo30', 'denovo31', 'denovo33', 'denovo38', 'denovo39', 'denovo51', 'denovo40', 'denovo42', 'denovo52', 'denovo45'])

        obs_hits_97 = pick_ref_contaminants(test_queries, test_ref_db_fp, test_input_fasta_fp, 0.97, test_output_dir)

        self.assertEqual(exp_hits_97, obs_hits_97)

        # Test at 99% ID level

        exp_hits_99 = set(['denovo20', 'denovo41', 'denovo40', 'denovo25', 'denovo24', 'denovo23', 'denovo46', 'denovo48', 'denovo0', 'denovo2', 'denovo4', 'denovo5', 'denovo6', 'denovo7', 'denovo12', 'denovo13', 'denovo10', 'denovo14', 'denovo15', 'denovo34', 'denovo36', 'denovo33', 'denovo52', 'denovo39', 'denovo51', 'denovo42'])

        obs_hits_99 = pick_ref_contaminants(test_queries, test_ref_db_fp, test_input_fasta_fp, 0.99, test_output_dir)

        self.assertEqual(exp_hits_99, obs_hits_99)

        # Test at 100% ID level

        exp_hits_100 = set(['denovo14', 'denovo25', 'denovo42', 'denovo20', 'denovo48'])

        obs_hits_100 = pick_ref_contaminants(test_queries, test_ref_db_fp, test_input_fasta_fp, 1.00, test_output_dir)

        self.assertEqual(exp_hits_100, obs_hits_100)

    def test_get_contamination_stats(self):
        """add_contamination_stats_to_biom: 
        """

        blank_sample_ids = ['Blank1', 'Blank2']

        test_biom = self.test_biom

        # test when passing already-proportional table

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, blank_sample_ids, proportional=True)

        # obs_contamination_header = exp_contamination_header
        # obs_contamination_stats_dict = exp_contamination_stats_dict

        # Header is as expected
        self.assertItemsEqual(exp_contamination_header,
                              obs_contamination_header)

        # Contamination stats dict is as expected
        self.assertDictEqual(exp_contamination_stats_dict,
                              obs_contamination_stats_dict)
        
        # test when passing non-proportional table

        (obs_contamination_header, obs_contamination_stats_dict) = get_contamination_stats(test_biom, blank_sample_ids, proportional=False)

        # obs_contamination_header = exp_contamination_header
        # obs_contamination_stats_dict = exp_contamination_stats_dict

        # Header is as expected
        self.assertItemsEqual(exp_contamination_header,
                              obs_contamination_header)

        # Contamination stats dict is as expected
        assertDeepAlmostEqual(self,exp_contamination_stats_dict_prop,
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

exp_contamination_header = ['maxS', 'avgS','maxB','avgB']
exp_contamination_stats_dict = {'otu1': [4,3,1,0.5], 'otu2': [6,3.5,4,2.5], 'otu3': [6,3,4,3.5], 'otu4': [3,1.5,0,0], 'contam1': [1,0.5,4,3], 'contam2': [3,1.5,4,2.5], 'contam3': [3,2.5,4,2], 'contam4': [0,0,3,2]}
exp_contamination_stats_dict_prop = {'otu1': [0.4,0.276923076923077,0.0416666666666667,0.0208333333333333], 'otu2': [0.230769230769231,0.215384615384615,0.166666666666667,0.145833333333333], 'otu3': [0.230769230769231,0.115384615384615,0.375,0.270833333333333], 'otu4': [0.115384615384615,0.0576923076923077,0,0], 'contam1': [0.0384615384615385,0.0192307692307692,0.25,0.208333333333333], 'contam2': [0.115384615384615,0.0576923076923077,0.166666666666667,0.145833333333333], 'contam3': [0.4,0.257692307692308,0.166666666666667,0.0833333333333333], 'contam4': [0,0,0.125,0.125]} 

test_sample_map = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tBlank\tAbund\tDescription\nBlank1\tAACCGCAT\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t1\t0\tBlank1\nBlank2\tAACCGCAC\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t1\t1\tBlank2\nSample1\tAACCGCAA\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t0\t2\tSample1\nSample2\tAACCGCAG\tCATGCTGCCTCCCGTAGGAGTGAGTTTGATCNTGGCTCAG\t0\t3\tSample2\n"""

# run tests if called from command line
if __name__ == "__main__":
    main()