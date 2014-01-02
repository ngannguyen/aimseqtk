#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.cdr3len.cdr3len
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.cdr3len.cdr3len as cdr3len
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone


class TestCdr3lenFuncs(unittest.TestCase):
    def setUp(self):
        c1 = lclone.Clone(1, 0.1, 'a', ['v1'], ['j1'])
        c1.__setitem__('cdr3aa', 'a')
        c2 = lclone.Clone(2, 0.2, 'abc', ['v1', 'v2'], ['j1'])
        c2.__setitem__('cdr3aa', 'abc')
        c3 = lclone.Clone(4, 0.4, 'bc', ['v3'], ['j1'])
        c3.__setitem__('cdr3aa', 'bc')
        c4 = lclone.Clone(3, 0.3, 'ab', ['v2'], ['j2'])
        c4.__setitem__('cdr3aa', 'ab')

        s1 = lsample.Sample('s1', [c1, c2, c3, c4])
        self.s1 = s1
        
        c4 = lclone.Clone(4, 0.2, 'abcdef', ['v1'], ['j1'])
        c4.__setitem__('cdr3aa', 'abcdef')
        c5 = lclone.Clone(16, 0.8, 'cde', ['v5'], ['j5'])
        c5.__setitem__('cdr3aa', 'cde')
        s2 = lsample.Sample('s2', [c4, c5])
        self.s2 = s2

    def test_sample_lendist_stat(self):
        stat = cdr3len.sample_lendist_stat(self.s1)
        l2c = {1: 0.25, 2: 0.5, 3: 0.25}
        l2r = {1: 0.1, 2: 0.7, 3: 0.2}
        self.assertEqual(stat.len2clones, l2c)
        self.assertEqual(stat.len2reads, l2r)
        self.assertEqual(stat.median_clones, 2)
        self.assertEqual(stat.median_reads, 2)

    def test_get_lens(self):
        stat1 = cdr3len.sample_lendist_stat(self.s1)
        stat2 = cdr3len.sample_lendist_stat(self.s2)
        ls = cdr3len.get_lens([stat1, stat2])
        self.assertEqual(ls, [1, 2, 3, 6])

if __name__ == '__main__':
    unittest.main()

