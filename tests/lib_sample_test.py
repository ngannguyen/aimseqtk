#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.lib.sample functions
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.clone as libclone


class TestFilterBySize(unittest.TestCase):
    '''Testing Sample filtering_by_size function
    '''
    def setUp(self):
        clone1 = libclone.Clone(1, 0.1, "ATCG", ["V1"], ["J1"])
        clone2 = libclone.Clone(5, 0.5, "CGTA", ["V2"], ["J2"])
        clone3 = libclone.Clone(2, 0.2, "AAAA", ["V3"], ["J3"])
        clone4 = libclone.Clone(2, 0.2, "TTTT", ["V4"], ["J4"])
        clones = [clone1, clone2, clone3, clone4]
        self.sample = libsample.Sample("Sample", clones, group="Group")

    def test_filtering(self):
        newsample = libsample.filter_by_size(self.sample, 0)
        self.assertEqual(newsample.numclone, 4)

        newsample = libsample.filter_by_size(self.sample, 2)
        self.assertEqual(newsample.numclone, 3)

        newsample = libsample.filter_by_size(self.sample, 2, 2)
        self.assertEqual(len(newsample.clones), 2)
        self.assertEqual(newsample.numclone, 2)
        self.assertEqual(newsample.size, 4)
        self.assertEqual(newsample.name, "Sample")
        self.assertEqual(newsample.group, "Group")
    
        newsample = libsample.filter_by_size(self.sample, 2, 2, 
                                             freqadjust=True)
        self.assertEqual(len(newsample.clones), 2)
        self.assertEqual(newsample.clones[0].freq, 2.0/4)

        newsample = libsample.filter_by_size(self.sample, 0.2, 0.2, True)
        self.assertEqual(len(newsample.clones), 2)

    def tearDown(self):
        self.sample = None

class TestSampling(unittest.TestCase):
    '''Testing Sample sampling function
    '''
    def setUp(self):
        clone1 = libclone.Clone(1, 0.1, "ATCG", ["V1"], ["J1"])
        clone2 = libclone.Clone(5, 0.5, "CGTA", ["V2"], ["J2"])
        clone3 = libclone.Clone(2, 0.2, "AAAA", ["V3"], ["J3"])
        clone4 = libclone.Clone(2, 0.2, "TTTT", ["V4"], ["J4"])
        clones = [clone1, clone2, clone3, clone4]
        self.sample = libsample.Sample("Sample", clones)
        self.s2v = {"ATCG": ["V1"], "CGTA": ["V2"], "AAAA": ["V3"], "TTTT": ["V4"]}
        self.s2j = {"ATCG": ["J1"], "CGTA": ["J2"], "AAAA": ["J3"], "TTTT": ["J4"]}

    def test_sampling(self):
        self.assertRaises(ValueError, libsample.sampling, None, 5)
        empty_sample = libsample.Sample("Empty", [])
        self.assertRaises(ValueError, libsample.sampling, empty_sample, 7)
        self.assertRaises(ValueError, libsample.sampling, self.sample, 0)
        self.assertRaises(ValueError, libsample.sampling, self.sample, 20)

        subsample = libsample.sampling(self.sample, 7)
        self.assertTrue(len(subsample.clones) <= 4)
        self.assertEqual(subsample.size, 7)
        size = sum([c.count for c in subsample.clones])
        self.assertEqual(size, 7)
        freq = sum([c.freq for c in subsample.clones])
        self.assertTrue(abs(1 - freq) < 0.01)
        
        seqs = [c.nuc for c in subsample.clones]
        uniqseqs = []
        for s in seqs:
            if s not in uniqseqs:
                uniqseqs.append(s)
        self.assertTrue(len(uniqseqs) == len(seqs))

        for c in subsample.clones:
            self.assertTrue(c.count > 0)
            self.assertTrue(c.freq > 0)
            self.assertTrue(c.nuc in self.s2v)
            self.assertTrue(c.vgenes == self.s2v[c.nuc])
            self.assertTrue(c.jgenes == self.s2j[c.nuc])

    def test_sampling_uniq(self):
        self.assertRaises(ValueError, libsample.sampling_uniq, None, 2)
        empty_sample = libsample.Sample("Empty", [])
        self.assertRaises(ValueError, libsample.sampling_uniq, empty_sample, 3)
        self.assertRaises(ValueError, libsample.sampling_uniq, self.sample, 0)
        self.assertRaises(ValueError, libsample.sampling_uniq, self.sample, 7)

        subsample = libsample.sampling_uniq(self.sample, 2)
        self.assertEqual(len(subsample.clones), 2)
        s2c = {"ATCG": 1, "CGTA": 5, "AAAA": 2, "TTTT": 2}
        s2f = {"ATCG": 0.1, "CGTA": 0.5, "AAAA": 0.2, "TTTT": 0.2}
        for c in subsample.clones:
            self.assertTrue(c.nuc in self.s2v)
            self.assertTrue(c.vgenes == self.s2v[c.nuc])
            self.assertTrue(c.jgenes == self.s2j[c.nuc])
            self.assertEqual(c.count, s2c[c.nuc])
            self.assertEqual(c.freq, s2f[c.nuc])


if __name__ == '__main__':
    unittest.main()

