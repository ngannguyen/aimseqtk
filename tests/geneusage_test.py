#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.geneusage
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.geneusage.geneusage as geneusage
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone


class TestGeneusageFuncs(unittest.TestCase):
    def setUp(self):
        c1 = lclone.Clone(1, 0.1, 'abc', ['v1'], ['j1'])
        c2 = lclone.Clone(5, 0.5, 'abc', ['v1', 'v2'], ['j1'])
        c3 = lclone.Clone(4, 0.4, 'bcd', ['v3'], ['j3'])
        s1 = lsample.Sample('s1', [c1, c2, c3])
        self.s1 = s1
        
        c4 = lclone.Clone(4, 0.2, 'abc', ['v1'], ['j1'], ['d1'])
        c5 = lclone.Clone(16, 0.8, 'cde', ['v5'], ['j5'])
        s2 = lsample.Sample('s2', [c4, c5])
        self.s2 = s2

    def test_sample_geneusage_stat(self):
        stat = geneusage.sample_geneusage_stat(self.s1)
        t2g2c = {'v': {'v1': 1.5/3, 'v2': 0.5/3, 'v3': 1.0/3},
                 'd': {},
                 'j': {'j1': 2.0/3, 'j3': 1.0/3},
                 'vj': {'v1|j1': 1.5/3, 'v2|j1': 0.5/3, 'v3|j3': 1.0/3},
                 'dj': {}}
        self.assertEqual(t2g2c, stat.type2gene2clones)
        
        t2g2r = {'v': {'v1': 0.35, 'v2': 0.25, 'v3': 0.4},
                 'd': {},
                 'j': {'j1': 0.6, 'j3': 0.4},
                 'vj': {'v1|j1': 0.35, 'v2|j1': 0.25, 'v3|j3': 0.4},
                 'dj': {}}
        self.assertEqual(t2g2r, stat.type2gene2reads)
        
        # s2 with d info
        stat2 = geneusage.sample_geneusage_stat(self.s2)
        t2g2c = {'v': {'v1': 1.0/2, 'v5': 1.0/2},
                 'd': {'d1': 1.0/2},
                 'j': {'j1': 1.0/2, 'j5': 1.0/2},
                 'vj': {'v1|j1': 1.0/2, 'v5|j5': 1.0/2},
                 'dj': {'d1|j1': 1.0/2}}
        self.assertEqual(t2g2c, stat2.type2gene2clones)
        
    def test_get_genes(self):
        stat1 = geneusage.sample_geneusage_stat(self.s1)
        stat2 = geneusage.sample_geneusage_stat(self.s2)
        vs = geneusage.get_genes([stat1, stat2], 'v')
        self.assertEqual(sorted(vs), ['v1', 'v2', 'v3', 'v5'])
        js = geneusage.get_genes([stat1, stat2], 'j')
        self.assertEqual(sorted(js), ['j1', 'j3', 'j5'])
        ds = geneusage.get_genes([stat1, stat2], 'd')
        self.assertEqual(sorted(ds), ['d1'])
        vjs = geneusage.get_genes([stat1, stat2], 'vj')
        self.assertEqual(sorted(vjs), ['v1|j1', 'v2|j1', 'v3|j3', 'v5|j5'])
        djs = geneusage.get_genes([stat1, stat2], 'dj')
        self.assertEqual(sorted(djs), ['d1|j1'])

    def test_get_geneusage(self):
        stat1 = geneusage.sample_geneusage_stat(self.s1)
        v = geneusage.get_geneusage(stat1, ('clones', 'v', 'v1'))
        self.assertEqual(v, 0.5)

        j = geneusage.get_geneusage(stat1, ('reads', 'vj', 'v2|j9'))
        self.assertEqual(j, 0.0)

if __name__ == '__main__':
    unittest.main()


