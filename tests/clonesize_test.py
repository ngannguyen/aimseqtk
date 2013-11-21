#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.clonesize
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.clonesize.repsize as repsize
import aimseqtk.lib.clone as libclone
import aimseqtk.lib.sample as libsample


class TestGettingGroupAvr(unittest.TestCase):
    def setUp(self):
        clone11 = libclone.Clone(9, 0.9, "ATCG", ["V1-1"], ["J1-1"])
        clone12 = libclone.Clone(1, 0.1, "CGTA", ["V1-2"], ["J1-2"])
        sample1 = libsample.Sample("S1", [clone11, clone12], "G1")
        
        clone31 = libclone.Clone(4, 0.2, "TTTT", ["V3-1"], ["J3-1"])
        clone32 = libclone.Clone(8, 0.4, "CCCC", ["V3-1"], ["J3-1"])
        clone33 = libclone.Clone(8, 0.4, "CCCC", ["V3-3"], ["J3-3"])
        sample21 = libsample.Sample("S21", [clone31, clone32, clone33], "G2")
        
        clone4 = libclone.Clone(100, 1.0, "GGG", ["V4"], ["J4"])
        sample22 = libsample.Sample("S22", [clone4], "G2")

        self.name2sample = {"S1": sample1, "S21": sample21, "S22": sample22}
        self.group2names = {"G1": ["S1"], "G2": ["S21", "S22"]}

    def test_get_group_avr(self):
        group2avr = repsize.get_group_avr(self.name2sample, self.group2names)
        self.assertEqual(len(group2avr), 2)
        groups = set(["G1", "G2"])
        self.assertEqual(set(group2avr.keys()), groups)
        
        avr1 = group2avr["G1"]
        size1 = 10
        numclone1 = 2
        self.assertEqual(avr1.name, "G1_Avr")
        self.assertEqual(avr1.group, "G1")
        self.assertEqual(avr1.size, size1)
        self.assertEqual(avr1.numclone, numclone1)

        avr2 = group2avr["G2"]
        size2 = 120/2
        numclone2 = (3 + 1)/2
        self.assertEqual(avr2.name, "G2_Avr")
        self.assertEqual(avr2.group, "G2")
        self.assertEqual(avr2.size, size2)
        self.assertEqual(avr2.numclone, numclone2)

class TestSortingSamplesByGroup(unittest.TestCase):
    def setUp(self):
        # Group 1
        clone11 = libclone.Clone(9, 0.9, "ATCG", ["V1-1"], ["J1-1"])
        clone12 = libclone.Clone(1, 0.1, "CGTA", ["V1-2"], ["J1-2"])
        sample1 = libsample.Sample("S1-1", [clone11, clone12], "G1")
        
        clone21 = libclone.Clone(50, 1.0, "AAAA", ["V2-1"], ["J2-1"])
        sample2 = libsample.Sample("S1-2", [clone21], "G1")

        avr1 = libsample.Sample("G1_Avr", [], "G1")
        avr1.size = 30
        avr1.numclone = 3/2

        # Group 2
        clone31 = libclone.Clone(4, 0.2, "TTTT", ["V3-1"], ["J3-1"])
        clone32 = libclone.Clone(8, 0.4, "CCCC", ["V3-1"], ["J3-1"])
        clone33 = libclone.Clone(8, 0.4, "CCCC", ["V3-3"], ["J3-3"])
        sample3 = libsample.Sample("S2-1", [clone31, clone32, clone33], "G2")
        
        avr2 = libsample.Sample("G2_Avr", [], "G2")
        avr2.size = 20
        avr2.numclone = 3

        # Group 3
        clone4 = libclone.Clone(100, 1.0, "GGG", ["V4"], ["J4"])
        sample4 = libsample.Sample("S3-1", [clone4], "G3")
        
        avr3 = libsample.Sample("G3_Avr", [], "G3")
        avr3.size = 100
        avr3.numclone = 1

        self.name2sample = {"S1-1": sample1, "S1-2": sample2,
                            "S2-1": sample3, "S3-1": sample4}
        self.group2names = {"G1": ["S1-1", "S1-2"], "G2": ["S2-1"],
                            "G3": ["S3-1"]}
        self.group2avr = {"G1": avr1, "G2": avr2, "G3": avr3}
        
    def test_sort_samples_by_group(self):
        samples = repsize.sort_samples_by_group(self.name2sample,
                                    self.group2names, True, self.group2avr)
        names = [s.name for s in samples]
        self.assertEqual(names, 
                         ["S1-2", "S1-1", "G1_Avr", "S2-1", "G2_Avr", "S3-1", "G3_Avr"])
        self.assertEqual(len(samples[1].clones), 2)
        self.assertEqual(samples[1].size, 10)
        self.assertEqual(samples[1].numclone, 2)
        self.assertEqual(samples[1].group, "G1")

        samples = repsize.sort_samples_by_group(self.name2sample,
                                    self.group2names, False, self.group2avr)
        names = [s.name for s in samples]
        self.assertEqual(names, ["S1-2", "S1-1", "S2-1", "S3-1"])


if __name__ == '__main__':
    unittest.main()

