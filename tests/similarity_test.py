#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.properties.similarity
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.properties.similarity as simi
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone


class TestGetRarefactionSizes(unittest.TestCase):
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

    def test_pair_similarity(self):
        attrs = ['numshare', 'chao', 'horn'] 
        stat = simi.pair_similarity(self.s1, self.s1, attrs, 'normcount')
        self.assertTrue(stat.numshare, 3)
        self.assertTrue(stat.chao, 1)
        self.assertTrue(stat.horn, 1)

        stat = simi.pair_similarity(self.s1, self.s2, attrs, 'normcount')
        self.assertTrue(stat.numshare == 1)
        self.assertTrue(stat.chao is not None)
        self.assertTrue(stat.horn is not None)

if __name__ == '__main__':
    unittest.main()

