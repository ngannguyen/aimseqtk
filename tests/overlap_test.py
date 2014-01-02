#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.overlap
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.overlap.overlap as overlap
import aimseqtk.src.overlap.trackclone as trackclone
import aimseqtk.lib.sample as lsample
import aimseqtk.lib.clone as lclone


class TestOverlapFuncs(unittest.TestCase):
    def setUp(self):
        self.names1 = ['s1', 's2', 's3']
        self.clone2sam2size = {'clone1': {'s1': 0.2, 's2': 0.1, 's3': 0.4, 's4': 0.1},
                          'clone2': {'s1': 0.1, 's3': 0.1, 's4': 0.1, 's5': 0.1, 's6': 0.1},
                          'clone3': {'s1': 0.1},
                          'clone4': {'s1': 0.3, 's2': 0.2},
                          'clone5': {'s4': 0.5, 's5': 0.2, 's6': 0.2}}

    def test_clone_group_freq(self):
        portion = overlap.clone_group_freq(self.clone2sam2size, self.names1, 'clone1')
        self.assertEqual(1.0, portion)
        portion = overlap.clone_group_freq(self.clone2sam2size, self.names1, 'clone2')
        self.assertEqual(2.0/3, portion)
        portion = overlap.clone_group_freq(self.clone2sam2size, self.names1, 'clone3')
        self.assertEqual(1.0/3, portion)
        portion = overlap.clone_group_freq(self.clone2sam2size, self.names1, 'clone5')
        self.assertEqual(0.0/3, portion)
    
    def test_group_major_clones(self):
        clones = overlap.group_major_clones(self.clone2sam2size, self.names1, 0.0)
        c = ['clone1', 'clone2', 'clone3', 'clone4', 'clone5']
        self.assertEqual(set(clones), set(c))
        
        clones = overlap.group_major_clones(self.clone2sam2size, self.names1, 0.5)
        c = ['clone1', 'clone2', 'clone4']
        self.assertEqual(set(clones), set(c))
        
        clones = overlap.group_major_clones(self.clone2sam2size, self.names1, 1.0)
        c = ['clone1']
        self.assertEqual(set(clones), set(c))

class TestTrackcloneFuncs(unittest.TestCase):
    def setUp(self):
        self.clone2sam2size = {'clone1': {'s11': 0.2, 's12': 0.1, 's21': 0.4, 's31': 0.1},
                          'clone2': {'s11': 0.1, 's31': 0.1, 's32': 0.1, 's21': 0.1, 's22': 0.1},
                          'clone3': {'s11': 0.1},
                          'clone4': {'s12': 0.3, 's22': 0.2},
                          'clone5': {'s22': 0.5, 's32': 0.2, 's12': 0.2}}

        self.groups = ['g1', 'g2', 'g3']
        self.group2samples = {'g1': ['s11', 's12'], 'g2': ['s21', 's22'], 'g3': ['s31', 's32']}
        self.group2samples0 = {'g1': ['s11', 's12'], 'g2': ['s21', 's22'], 'g3': ['s31', 's32'], 'g': ['s0']}
        self.sam2group = {'s11': 'g1', 's12': 'g1', 's21': 'g2', 's22': 'g2', 's31': 'g3', 's32': 'g3', 's0': 'g'}

    def test_get_matched_samples(self):
        sams, i = trackclone.get_matched_samples('s11', 'g1', self.groups, self.group2samples)
        self.assertEqual(sams, ['s11', 's21', 's31'])
        self.assertEqual(i, 0)
        sams, i = trackclone.get_matched_samples('s22', 'g2', self.groups, self.group2samples)
        self.assertEqual(sams, ['s12', 's22', 's32'])
        self.assertEqual(i, 1)
        sams, i = trackclone.get_matched_samples('s31', 'g3', ['g1', 'g2'], self.group2samples)
        self.assertEqual(sams, ['s11', 's21'])

        self.assertRaises(ValueError, trackclone.get_matched_samples, 's22', 'g1', self.groups, self.group2samples)
        self.assertRaises(ValueError, trackclone.get_matched_samples, 's31', 'g3', ['g1', 'g4'], self.group2samples)
        self.assertRaises(ValueError, trackclone.get_matched_samples, 's32', 'g3', ['g1', 'g'], self.group2samples0)
        
    def test_track_clone(self):
        rows, indices = trackclone.track_clone('clone1', self.clone2sam2size['clone1'], self.groups, self.group2samples, self.sam2group)
        e_rows = [[0.2, 0.4, 0.1], [0.1, 0.0, 0.0]]
        e_indices = [0, 1]
        self.assertEqual(rows, e_rows)
        self.assertEqual(indices, e_indices)
    
class TestTopclonesFuncs(unittest.TestCase):
    def setUp(self):
        c1 = lclone.Clone(1, 0.1, 'a', ['v1'], ['j1'])
        c2 = lclone.Clone(2, 0.2, 'abc', ['v1', 'v2'], ['j1'])
        c3 = lclone.Clone(4, 0.4, 'bc', ['v3'], ['j1'])
        c4 = lclone.Clone(3, 0.3, 'ab', ['v2'], ['j2'])

        s1 = lsample.Sample('s1', [c1, c2, c3, c4])
        self.s1 = s1
        
        c4 = lclone.Clone(4, 0.2, 'abc', ['v1'], ['j1'])
        c5 = lclone.Clone(8, 0.4, 'cde', ['v5'], ['j5'])
        c6 = lclone.Clone(8, 0.4, 'ab', ['v2'], ['j2'])
        s2 = lsample.Sample('s2', [c4, c5, c6])
        self.s2 = s2
        
        self.sams = [self.s1, self.s2]
        self.groups = ['g1', 'g2']
        self.g2s = {'g1': ['s1'], 'g2': ['s2']}
        self.s2g = {'s1': 'g1', 's2': 'g2'}

    def test_sample_top_clones(self):
        topclone2size = trackclone.sample_top_clones(self.s1, 0.0, attr='freq')
        self.assertEqual(len(topclone2size), 5)
        topclone2size = trackclone.sample_top_clones(self.s1, 2, attr='count')
        tc2s = {'v3_bc_j1': 4, 'v2_ab_j2': 3}
        self.assertEqual(topclone2size, tc2s)
        topclone2size = trackclone.sample_top_clones(self.s1, 0.15)
        tc2s = {'v3_bc_j1': 0.4, 'v2_ab_j2': 0.3}
        self.assertEqual(topclone2size, tc2s)
    
    def test_top_clones(self):
        tc2n2s = trackclone.top_clones(self.sams, 0.2) 
        e_tc2n2s = {'v3_bc_j1': {'s1': 0.4}, 'v2_ab_j2': {'s1': 0.3, 's2': 0.4},
                    'v1_abc_j1': {'s2': 0.2}, 'v5_cde_j5': {'s2': 0.4}}
        self.assertEqual(tc2n2s, e_tc2n2s)
        
    def test_clone_get_samples(self):
        n2s = trackclone.clone_get_samples('v1_abc_j1', self.sams, sizetype='freq')
        e_n2s = {'s1': 0.1, 's2': 0.2}
        self.assertEqual(n2s, e_n2s)

        n2s = trackclone.clone_get_samples('v1_abc_j7', self.sams, sizetype='freq')
        self.assertEqual(n2s, {})

    def test_track_top_clones(self):
        tc2rows, tc2indices = trackclone.track_top_clones(self.sams, 0.2, self.groups, self.g2s, self.s2g)
        e_tc2rows = {'v3_bc_j1': [[0.4, 0.0]], 'v2_ab_j2': [[0.3, 0.4]],
                     'v1_abc_j1': [[0.1, 0.2]], 'v5_cde_j5': [[0.0, 0.4]]}
        e_tc2indices = {'v3_bc_j1': [0], 'v2_ab_j2': [0],
                     'v1_abc_j1': [0], 'v5_cde_j5': [0]}
        self.assertEqual(tc2rows, e_tc2rows)
        self.assertEqual(tc2indices, e_tc2indices)

if __name__ == '__main__':
    unittest.main()


