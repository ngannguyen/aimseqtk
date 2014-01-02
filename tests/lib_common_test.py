#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.lib.common functions
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.lib.common as lcommon
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone

def myfunc(item):
    sample = item[0]
    return sample.numclone

class TestCommonFunctions(unittest.TestCase):
    def setUp(self):
        dict1 = {'a': 5, 'b': -1, 'c': 2, 'd': 10}
        self.dict1 = dict1

        sam1 = lsample.Sample('sam1', group='g1')
        sam1.size = 10
        sam1.numclone = 2
        self.sam1 = sam1
        
        sam2 = lsample.Sample('sam2', group='g2')
        sam2.size = 5
        sam2.numclone = 4
        self.sam2 = sam2

        sam3 = lsample.Sample('sam3', group='g1')
        sam3.size = 15
        sam3.numclone = 1
        self.sam3 = sam3

        dict2 = {'sam1': sam1, 'sam2': sam2, 'sam3': sam3}
        self.dict2 = dict2
         
        dict3 = {'sam12': [1, 2], 'sam23': [2, 3]}
        self.dict3 = dict3
         
    def test_get_cumulative(self):
        vals = [1, 0, -1, 0.5, 2.5]
        forward_vals = [1, 1, 0, 0.5, 3]
        cumu_vals = lcommon.get_cumulative(vals, forward=True)
        self.assertEqual(forward_vals, cumu_vals)
        backward_vals = [3, 2, 2, 3, 2.5] 
        cumu_vals = lcommon.get_cumulative(vals)
        self.assertEqual(backward_vals, cumu_vals)
        
    def test_sort_dict_by_value(self):
        # test normal dictionary
        sorted_dict1 = lcommon.sort_dict_by_value(self.dict1)
        dict1 = [('d', 10), ('a', 5), ('c', 2), ('b', -1)]
        self.assertEqual(sorted_dict1, dict1)

        sorted_dict2 = lcommon.sort_dict_by_value(self.dict2)
        dict2 = [('sam3', self.sam3), ('sam1', self.sam1), ('sam2', self.sam2)]
        self.assertEqual(sorted_dict2, dict2)

        sorted_dict2 = lcommon.sort_dict_by_value(self.dict2, myfunc)
        dict2 = [('sam2', self.sam2), ('sam1', self.sam1), ('sam3', self.sam3)]
        self.assertEqual(sorted_dict2, dict2)

    def test_get_val2keys(self):
        v2ks = lcommon.get_val2keys(self.dict3)
        for v, ks in v2ks.iteritems():
            v2ks[v] = sorted(ks)
        dict3 = {1: ['sam12'], 2: ['sam12', 'sam23'], 3: ['sam23']}
        self.assertEqual(v2ks, dict3)
        self.assertRaises(lcommon.InputError, lcommon.get_val2key_1to1, self.dict3)

class TestGroupCommonFunctions(unittest.TestCase):
    def setUp(self):
        sam1 = lsample.Sample('sam1', group='g1')
        sam1.size = 10
        sam1.numclone = 2
        sam1.__setitem__('clones', [5, 0])
        self.sam1 = sam1
        
        sam2 = lsample.Sample('sam2', group='g1')
        sam2.size = 5
        sam2.numclone = 4
        sam2.__setitem__('clones', [3, 2])
        self.sam2 = sam2

        sam3 = lsample.Sample('sam3', group='g2')
        sam3.size = 15
        sam3.numclone = 1
        sam3.__setitem__('clones', ['a', 'b', 'c'])
        self.sam3 = sam3

        name2obj = {'sam1': sam1, 'sam2': sam2, 'sam3': sam3}
        self.name2obj = name2obj

        group2names = {'g1': ['sam1', 'sam2'], 'g2': ['sam3']}
        self.group2names = group2names

    def test_get_group_avr(self):
        group2avr = lcommon.get_group_avr(self.name2obj, self.group2names)
        self.assertEqual(len(group2avr), 2)
        self.assertEqual(sorted(group2avr.keys()), ['g1', 'g2'])
        (name1, avr1) = group2avr['g1']
        self.assertEqual(name1, 'g1_Avr')
        self.assertEqual(avr1.size, (10 + 5) / 2)
        self.assertEqual(avr1.numclone, (2 + 4) / 2)
        self.assertEqual(avr1.clones, [(5 + 3) / 2, (0 + 2) / 2])

        (name2, avr2) = group2avr['g2']
        self.assertEqual(name2, 'g2_Avr')
        self.assertEqual(avr2.size, 15)
        self.assertEqual(avr2.numclone, 1)
        self.assertEqual(avr2.clones, None)

    def test_sort_objs_by_group(self):
        group2avr = lcommon.get_group_avr(self.name2obj, self.group2names)
        sorted_objs = lcommon.sort_objs_by_group(self.name2obj,
                                        self.group2names, True, group2avr)
        self.assertEqual(len(sorted_objs), 5)
        names = ['sam1', 'sam2', 'g1_Avr', 'sam3', 'g2_Avr']
        self.assertEqual(names, [o[0] for o in sorted_objs])
        n2o = self.name2obj
        objs = [n2o['sam1'], n2o['sam2'], group2avr['g1'][1], n2o['sam3'],
                                                        group2avr['g2'][1]]
        self.assertEqual(objs, [o[1] for o in sorted_objs])
        #Without adding average group
        sorted_objs = lcommon.sort_objs_by_group(self.name2obj,
                                        self.group2names, False, group2avr)
        self.assertEqual(len(sorted_objs), 3)
        names = ['sam1', 'sam2', 'sam3']
        self.assertEqual(names, [o[0] for o in sorted_objs])
        n2o = self.name2obj
        objs = [n2o['sam1'], n2o['sam2'], n2o['sam3']]
        self.assertEqual(objs, [o[1] for o in sorted_objs])

if __name__ == '__main__':
    unittest.main()

