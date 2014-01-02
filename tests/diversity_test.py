#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.properties
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.properties.diversity as diversity
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.clone as libclone
import aimseqtk.lib.common as libcommon


class TestGetRarefactionSizes(unittest.TestCase):
    def setUp(self):
        sample = libsample.Sample('s')
        sample.size = 100
        sample.numclone = 11
        self.sample = sample

    def test_sample_rf_sizes(self):
        sizes = diversity.sample_rf_sizes(self.sample, 20)
        expected = [20, 40, 60, 80, 100]
        self.assertEqual(set(expected), set(sizes))

        sizes = diversity.sample_rf_sizes(self.sample, 20, [10, 80, 150])
        expected = [10, 80]
        self.assertEqual(set(expected), set(sizes))

        sizes = diversity.sample_rf_sizes(self.sample, rf_sizes=[200, 500])
        expected = [100]
        self.assertEqual(set(expected), set(sizes))
        
        sizes = diversity.sample_rf_sizes(self.sample)
        expected = [100]
        self.assertEqual(set(expected), set(sizes))

class TestSamplingDiversity(unittest.TestCase):
    def setUp(self):
        clone1 = libclone.Clone(1, 0.1, "ATCG", ["V1"], ["J1"])
        clone2 = libclone.Clone(5, 0.5, "CGTA", ["V2"], ["J2"])
        clone3 = libclone.Clone(2, 0.2, "AAAA", ["V3"], ["J3"])
        clone4 = libclone.Clone(2, 0.2, "TTTT", ["V4"], ["J4"])
        clones = [clone1, clone2, clone3, clone4]
        self.sample = libsample.Sample("Sample", clones)
        
    def test_sampling_diversity(self):
        indices = ['numclone', 'fisher_alpha', 'simpson', 'shannon',
                   'invsimpson']
        sampling = diversity.sample_sampling_diversity(self.sample, 
                                                       [5, indices])
        for i in indices:
            self.assertTrue(i in sampling.getitems())
            self.assertTrue(sampling[i] is not None)
        numclone = sampling['numclone']
        self.assertTrue(numclone <= 4)

if __name__ == '__main__':
    unittest.main()


