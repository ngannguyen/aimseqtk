#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''Test aimseqtk.normalize.normalize functions
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.normalize.normalize as norm
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone
import aimseqtk.lib.common as lcommon


class TestMetagenomeSeqNormalization(unittest.TestCase):
    def setUp(self):
        c1 = lclone.Clone(1, 0.1, 'abc', ['v1'], ['j1'])
        c2 = lclone.Clone(5, 0.5, 'abc', ['v1', 'v2'], ['j1'])
        c3 = lclone.Clone(4, 0.4, 'bcd', ['v3'], ['j1'])
        s1 = lsample.Sample('s1', [c1, c2, c3])
        self.s1 = s1
        
        c4 = lclone.Clone(4, 0.2, 'abc', ['v1'], ['j1'])
        c5 = lclone.Clone(16, 0.8, 'cde', ['v5'], ['j5'])
        s2 = lsample.Sample('s2', [c4, c5])
        self.s2 = s2

        self.samples = [s1, s2]
        self.name2obj = {'s1': s1, 's2': s2}

    def test_clone_matrix(self):
        # cols = samples, rows = clones
        rows, colnames, rownames = norm.clone_matrix(self.samples, 'normcount')
        e_colnames = ['s1', 's2']
        self.assertEqual(e_colnames, colnames)
        e_rownames = ['v1_abc_j1', 'v2_abc_j1', 'v3_bcd_j1', 'v5_cde_j5']
        self.assertEqual(set(e_rownames), set(rownames))
        e_rows = [[3, 4], [2, 0], [4, 0], [0, 16]]
        e_n2r = {}
        for i, n in enumerate(e_rownames):
            e_n2r[n] = e_rows[i]
        e_rows = []
        for n in rownames:
            e_rows.extend(e_n2r[n])
        self.assertEqual(rows, e_rows)

    def test_metagenome_seq_norm(self):
        norm_matrix = norm.normalize_MRexp(self.samples, "normcount")
        self.assertEqual(len(norm_matrix), 8)
        samples = norm.matrix_to_normcount(norm_matrix, self.samples)
        names = ['s1', 's2']
        self.assertEqual(names, [s.name for s in samples])
        s1 = samples[0]
        c1 = s1.clones[0]
        self.assertTrue(abs(c1.normcount - 1500) / 1500 < 0.001)
        self.assertTrue(abs(c1['normcount'] -1500) / 1500 < 0.001)
        c2 = s1.clones[1]
        self.assertTrue(abs(c2.normcount -2500) / 2500 < 0.001)

if __name__ == '__main__':
    unittest.main()
