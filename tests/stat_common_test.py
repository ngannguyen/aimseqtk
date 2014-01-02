#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.lib.statcommon functions
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.lib.statcommon as stat
import aimseqtk.lib.common as lcommon
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone


def myfunc(obj, args):
    val = 0
    for attr in args:
        val += obj[attr]
    return val

class TestStatCommonFuncs(unittest.TestCase):
    def setUp(self):
        c1 = lclone.Clone(1, 0.1, 'abc', ['v1'], ['j1'])
        c2 = lclone.Clone(5, 0.5, 'abc', ['v1', 'v2'], ['j1' ,'j2'])
        c3 = lclone.Clone(4, 0.4, 'bcd', ['v3'], ['j1'])
        s1 = lsample.Sample('s1', [c1, c2, c3])
        self.s1 = s1
        
        c4 = lclone.Clone(4, 0.2, 'abc', ['v1'], ['j1'])
        c5 = lclone.Clone(16, 0.8, 'cde', ['v5'], ['j5'])
        s2 = lsample.Sample('s2', [c4, c5])
        self.s2 = s2

        self.samples = [s1, s2]
        self.name2obj = {'s1': s1, 's2': s2}

    def test_get_clone2samples(self):
        c2s2size = stat.get_clone2samples(self.samples)
        e_c2s2size = {'v1_abc_j1': {'s1': 2, 's2': 4}, 
                      'v1_abc_j2': {'s1': 1},
                      'v2_abc_j1': {'s1': 1},
                      'v2_abc_j2': {'s1': 1},
                      'v3_bcd_j1': {'s1': 4},
                      'v5_cde_j5': {'s2': 16}
                     }
        self.assertEqual(e_c2s2size, c2s2size)
        
        c2s2size = stat.get_clone2samples(self.samples, 'normcount')
        e_c2s2size = {'v1_abc_j1': {'s1': 2, 's2': 4}, 
                      'v1_abc_j2': {'s1': 1},
                      'v2_abc_j1': {'s1': 1},
                      'v2_abc_j2': {'s1': 1},
                      'v3_bcd_j1': {'s1': 4},
                      'v5_cde_j5': {'s2': 16}
                     }
        self.assertEqual(e_c2s2size, c2s2size)
        
        c2s2size = stat.get_clone2samples(self.samples, 'freq')
        e_c2s2size = {'v1_abc_j1': {'s1': 0.1 + 0.5/4, 's2': 0.2}, 
                      'v1_abc_j2': {'s1': 0.5/4},
                      'v2_abc_j1': {'s1': 0.5/4},
                      'v2_abc_j2': {'s1': 0.5/4},
                      'v3_bcd_j1': {'s1': 0.4},
                      'v5_cde_j5': {'s2': 0.8}
                     }
        self.assertEqual(e_c2s2size, c2s2size)
        
        self.assertRaises(KeyError, stat.get_clone2samples, self.samples, 'a')
        
    def test_group_vector(self):
        names = ['s1', 's2']
        v = stat.group_vector(names, self.name2obj)
        self.assertEqual([self.s1, self.s2], v)

        v = stat.group_vector(names, self.name2obj, 'size')
        self.assertEqual([10, 20], v)
        
        v = stat.group_vector(names, self.name2obj, None, myfunc,
                                                        ('size', 'numclone'))
        self.assertEqual([13, 22], v)

if __name__ == '__main__':
    unittest.main()


