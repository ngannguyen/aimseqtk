#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.recomb.aa_events
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.recomb.aa_events as aaevents
import aimseqtk.lib.clone as lclone

class TestAaeventsFuncs(unittest.TestCase):
    #def setUp(self):
     
    def test_get_all_nts(self):
        codon_lists = [['a', 'b'], ['c'], ['d', 'efg']]
        expected = ['acd', 'acefg', 'bcd', 'bcefg']
        received = aaevents.get_all_nts(codon_lists)
        self.assertEqual(set(expected), set(received))
        #self.assertRaises
    
    def test_left_max_match(self):
        seq1 = 'ATCGGG'
        seq2 = 'ATCGCCATCGGT'
        expected = 'ATCG'
        received = aaevents.left_max_match(seq1, seq2)
        self.assertEqual(expected, received)

    def test_right_max_match(self):
        seq1 =      'ATCGTGCA'
        seq2 = 'ATCGTTGCCCGCA'
        e = 'GCA'
        r = aaevents.right_max_match(seq1, seq2)
        self.assertEqual(e, r)
        seq3 = 'AAAAAAATCGTATCGTGCA'
        r3 = aaevents.right_max_match(seq1, seq3)
        self.assertEqual(seq1, r3)

    def test_find_min_vdel(self):
        vnt = 'tgtgccagcagcttaga'.upper()
        cdr3aa = 'CASSLTGCHQTF'
        e = 2
        r = aaevents.find_min_vdel(vnt, cdr3aa)
        self.assertEqual(e, r)
        cdr3aa = 'CASSLDATG'
        r = aaevents.find_min_vdel(vnt, cdr3aa)
        self.assertEqual(0, r)
        cdr3aa = 'CASSLEATG'
        r = aaevents.find_min_vdel(vnt, cdr3aa)
        self.assertEqual(0, r)
        cdr3aa = 'CASSSVATG'
        r = aaevents.find_min_vdel(vnt, cdr3aa)
        self.assertEqual(4, r)

    def test_find_min_jdel(self):
        jnt = 'tgaacactgaagctttcttt'.upper()
        cdr3aa = 'NTEAFF'
        r = aaevents.find_min_jdel(jnt, cdr3aa)
        self.assertEqual(2, r)

        cdr3aa = 'HKAFF'
        r = aaevents.find_min_jdel(jnt, cdr3aa)
        self.assertEqual(9, r)

        cdr3aa = 'HVAFF'
        r = aaevents.find_min_jdel(jnt, cdr3aa)
        self.assertEqual(10, r)

    def test_find_dmatches(self):
        daa = 'GCT'
        cdr3aa = 'CASSHTLGCTF'
        matches = aaevents.find_dmatches(daa, cdr3aa)
        e = [(7, 10)]
        self.assertEqual(e, matches)

        daa = 'AGCT'
        cdr3aa = 'CASSHTLGCTF'
        matches = aaevents.find_dmatches(daa, cdr3aa)
        e = []
        self.assertEqual(e, matches)

        daa = 'GC'
        cdr3aa = 'CASSGCHTLGCTF'
        matches = aaevents.find_dmatches(daa, cdr3aa)
        e = [(4, 6), (9, 11)]
        self.assertEqual(set(e), set(matches))
    
    def test_find_devents(self):
        d_nt = 'GGGACA'   # frame0 = GT, frame1 = G, frame2 = D
        cdr3aa = 'AGTT'
        
        devents = aaevents.find_devents(d_nt, cdr3aa)
        #for d in devents:
        #    print "%d\t%d\t*%s*\t%d\t%d\t*%s*" % (d.d5del, d.d3del, d.left_nts, d.cdr3aa_dstart, d.cdr3aa_dend, d.right_nts)
        #print len(devents)
        #d1 = aaevents.Devent(0, 0, '', 1, 3, '')
    
    def test_get_vdins_events(self):
        vnt = 'tgtgccagcagcttaga'.upper()
        cdr3aa = 'CASSLAGTQTF'
        vdel = 2
        devent = aaevents.Devent(0, 0, '', 6, 8, '')
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        e = ['AGCA', 'AGCC', 'AGCG', 'AGCT']
        self.assertEqual(set(e), set(vdins))

        vdel = 0
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        self.assertTrue(vdins is None)
        
        cdr3aa = 'CASSLDGTQTF'
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        e = ['AT', 'AC']
        self.assertEqual(set(e), set(vdins))

        cdr3aa = 'CASSLDAGTQTF'
        devent = aaevents.Devent(0, 0, '', 7, 9, '')
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        e = ['ATGCT', 'ATGCC', 'ATGCA', 'ATGCG', 'ACGCT', 'ACGCC', 'ACGCA', 'ACGCG']
        self.assertEqual(set(e), set(vdins))

        cdr3aa = 'CASSLDGTQTF'
        devent = aaevents.Devent(1, 0, 'T', 6, 8, '')
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        e = ['A']
        self.assertEqual(set(e), set(vdins))

        cdr3aa = 'CASSLDGTQTF'
        devent = aaevents.Devent(1, 0, 'TGG', -1, -1, '')
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        self.assertTrue( vdins is None)

        cdr3aa = 'CASSLDAGTQTF'
        devent = aaevents.Devent(1, 3, 'A', 7, 9, '')
        vdins = aaevents.get_vdins_events(vdel, vnt, devent, cdr3aa)
        e = ['ATGC', 'ACGC']
        self.assertEqual(set(e), set(vdins))

    def test_get_djins_events(self):
        #jnt = 'tga ac act gaa gct ttc ttt'.upper()
        jnt = 'tgaacactgaagctttcttt'.upper()
        cdr3aa = 'AGTVNTEAFF'
        devent = aaevents.Devent(2, 4, '', 1, 3, '')
        jdel = 0
        djins = aaevents.get_djins_events(jdel, jnt, devent, cdr3aa)
        e = ['GT']
        self.assertEqual(set(e), set(djins))

        #cdr3aa = 'AGTVN TEAFF'
        devent = aaevents.Devent(2, 4, '', 1, 3, '')
        jdel = 3
        djins = aaevents.get_djins_events(jdel, jnt, devent, cdr3aa)
        e = ['GTTAA', 'GTCAA', 'GTAAA', 'GTGAA']
        self.assertEqual(set(e), set(djins))

        devent = aaevents.Devent(2, 4, '', 1, 3, 'GT')
        jdel = 3
        djins = aaevents.get_djins_events(jdel, jnt, devent, cdr3aa)
        e = ['TAA', 'CAA', 'AAA', 'GAA']
        self.assertEqual(set(e), set(djins))

    
    def test_get_vjins_emptyd(self):
        vnt = 'tgtgccagcagcttaga'.upper()
        jnt = 'tgaacactgaagctttcttt'.upper()
        cdr3aa = 'CASSLDVNTEAFF'
        vdel = 0
        jdel = 0
        dnts = '' 
        e = ['ATGT', 'ACGT']
        vjins = aaevents.get_vjins_emptyd(vnt, vdel, jnt, jdel, dnts, cdr3aa)
        self.assertEqual(set(e), set(vjins))

        cdr3aa = 'CASSLDNVNTEAFF'
        e = ['ATAATGT', 'ACAATGT', 'ATAACGT', 'ACAACGT']
        vjins = aaevents.get_vjins_emptyd(vnt, vdel, jnt, jdel, dnts, cdr3aa)
        self.assertEqual(set(e), set(vjins))
        
        dnts = 'TAA'
        e = ['ATAATGT', 'ATAACGT']
        vjins = aaevents.get_vjins_emptyd(vnt, vdel, jnt, jdel, dnts, cdr3aa)
        self.assertEqual(set(e), set(vjins))

if __name__ == '__main__':
    unittest.main()


