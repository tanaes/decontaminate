#!/usr/bin/env python

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"

from unittest import TestCase, main

import numpy
from numpy import array
from numpy.testing import assert_almost_equal, assert_allclose

from biom.parse import parse_biom_table

from qiime.decontaminate import get_contamination_stats, compare_blank_abundances

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

        test_biom = parse_biom_table(open(test_biom_fp,'Ur'))

    def test_get_contamination_stats(self):
        """add_contamination_stats_to_biom: 
        """
        test_biom_fp = '/Users/jonsanders/Development/git_sw/qiime/qiime_test_data/decontaminate/test_otu_table.biom'

        test_biom = parse_biom_table(open(test_biom_fp,'Ur'))
                # print('testing contamination')

        blank_sample_ids = ['Blank1', 'Blank2']

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

        # test for maxS > maxB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 exp_contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'maxB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu3','otu4'],passed_seqs)


        # test for maxS > avgB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 exp_contamination_header,
                                 sample_stat = 'maxS',
                                 blank_stat = 'avgB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu3','otu4','contam2','contam3'],passed_seqs)


        # test for avgS > avgB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 exp_contamination_header,
                                 sample_stat = 'avgS',
                                 blank_stat = 'avgB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu2','otu4','contam3'],passed_seqs)


        # test for avgS > maxB
        passed_seqs = compare_blank_abundances(contamination_stats_dict,
                                 exp_contamination_header,
                                 sample_stat = 'avgS',
                                 blank_stat = 'maxB',
                                 scalar = 1,
                                 negate = False)
        
        self.assertItemsEqual(['otu1','otu4'],passed_seqs)


test_biom_file = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","matrix_type": "sparse","generated_by": "BIOM-Format 2.0.1-dev","date": "2014-09-14T17:26:13.141629","type": "OTU table","matrix_element_type": "float","shape": [8, 4],"data": [[0,0,1.0],[0,2,4.0],[0,3,2.0],[1,0,4.0],[1,1,1.0],[1,2,6.0],[1,3,1.0],[2,0,4.0],[2,1,3.0],[2,2,6.0],[3,2,3.0],[4,0,4.0],[4,1,2.0],[4,2,1.0],[5,0,4.0],[5,1,1.0],[5,2,3.0],[6,0,4.0],[6,2,3.0],[6,3,2.0],[7,0,3.0],[7,1,1.0]],"rows": [{"id": "otu1", "metadata": null},{"id": "otu2", "metadata": null},{"id": "otu3", "metadata": null},{"id": "otu4", "metadata": null},{"id": "contam1", "metadata": null},{"id": "contam2", "metadata": null},{"id": "contam3", "metadata": null},{"id": "contam4", "metadata": null}],"columns": [{"id": "Blank1", "metadata": null},{"id": "Blank2", "metadata": null},{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null}]}"""

exp_contamination_header = ['maxS', 'avgS','maxB','avgB']
exp_contamination_stats_dict = {'otu1': [4,3,1,0.5], 'otu2': [6,3.5,4,2.5], 'otu3': [6,3,4,3.5], 'otu4': [3,1.5,0,0], 'contam1': [1,0.5,4,3], 'contam2': [3,1.5,4,2.5], 'contam3': [3,2.5,4,2], 'contam4': [0,0,3,2]}
exp_contamination_stats_dict_prop = {'otu1': [0.4,0.276923076923077,0.0416666666666667,0.0208333333333333], 'otu2': [0.230769230769231,0.215384615384615,0.166666666666667,0.145833333333333], 'otu3': [0.230769230769231,0.115384615384615,0.375,0.270833333333333], 'otu4': [0.115384615384615,0.0576923076923077,0,0], 'contam1': [0.0384615384615385,0.0192307692307692,0.25,0.208333333333333], 'contam2': [0.115384615384615,0.0576923076923077,0.166666666666667,0.145833333333333], 'contam3': [0.4,0.257692307692308,0.166666666666667,0.0833333333333333], 'contam4': [0,0,0.125,0.125]} 


# run tests if called from command line
if __name__ == "__main__":
    main()