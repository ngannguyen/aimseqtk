#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Testing aimseqtk.src.input functions
'''

import os
import sys
import unittest2 as unittest

import aimseqtk.src.input.inputcommon as inputcommon
import aimseqtk.src.input.sequenta as sequenta
import aimseqtk.src.input.adaptive as adaptive
import aimseqtk.src.input.mitcr as mitcr


class TestMitcrParseline(unittest.TestCase):
    '''Testing MiTCR parsing function
    '''
    def setUp(self):
        self.index2col = {0: "Read count", 1: "Percentage", 
                2: "CDR3 nucleotide sequence", 3: "CDR3 nucleotide quality",
                4: "Min quality", 5: "CDR3 amino acid sequence",
                6: "V alleles",
                7: "V segments", 8: "J alleles", 9: "J segments",
                10: "D alleles",
                11: "D segments", 12: "Last V nucleotide position",
                13: "First D nucleotide position",
                14: "Last D nucleotide position",
                15: "First J nucleotide position", 16: "VD insertions",
                17: "DJ insertions", 18: "Total insertions"}
        self.productive = ("11763\t0.13098233971004164\tTGTGCCAGCAGCTTAGGGGAA"
                           + "AACATTCAGTACTTC\tJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"
                           + "JJJJJJ\t41\tCASSLGENIQYF\tTRBV13*01\tTRBV13\t"
                           + "TRBJ2-4*01\tTRBJ2-4\t"
                           + "TRBD1*01, TRBD2*02, TRBD2*01\tTRBD1, TRBD2\t16"
                           + "\t-1\t-1\t19\t-1\t-1\t2")

    def test_parse_productive_clone(self):
        clone = mitcr.mitcr_parseline(self.productive, self.index2col)
        self.assertTrue(clone is not None)
        self.assertEqual(clone.count, 11763)
        self.assertEqual("%.3f" % (clone.freq * 100), "0.131")
        
        vgenes = ["TRBV13"]
        self.assertTrue(set(clone.vgenes) == set(vgenes))
        valleles = ["TRBV13*01"]
        self.assertTrue(set(clone.valleles) == set(valleles))
        jgenes = ["TRBJ2-4"]
        self.assertTrue(set(clone.jgenes) == set(jgenes))
        jalleles = ["TRBJ2-4*01"]
        self.assertTrue(set(clone.jalleles) == set(jalleles))
        dgenes = ["TRBD1", "TRBD2"]
        self.assertTrue(set(clone.dgenes) == set(dgenes))
        dalleles = ["TRBD1*01", "TRBD2*01", "TRBD2*02"]
        self.assertTrue(set(clone.dalleles) == set(dalleles))

        self.assertTrue(clone.productive)
        nuc = ("TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC")
        self.assertEqual(clone.nuc, nuc)
        self.assertEqual(clone.cdr3nuc, nuc)
        aa = "CASSLGENIQYF"
        self.assertEqual(clone.aa, aa)
        self.assertEqual(clone.cdr3aa, aa)

        self.assertEqual(clone.lastvpos, 16)
        self.assertEqual(clone.firstdpos,  -1)
        self.assertEqual(clone.lastdpos, -1)
        self.assertEqual(clone.firstjpos, 19)
         
class TestAdaptiveParseline(unittest.TestCase):
    '''Testing Adaptive parsing function
    '''
    def setUp(self):
        self.index2col = {0: 'sequenceID', 1: 'container', 2: 'nucleotide',
                           3: 'aminoAcid', 4: 'normalizedFrequency', 
                           5: 'normalizedCopy', 6: 'rawFrequency', 7: 'copy',
                           8: 'cdr3Length', 9: 'VFamilyName', 10: 'VGeneName',
                           11: 'VTies', 12: 'DGeneName', 13: 'JGeneName',
                           14: 'JTies', 15: 'VDeletion', 16: 'd5Deletion',
                           17: 'd3Deletion', 18: 'JDeletion',
                           19: 'n2Insertion', 20: 'n1Insertion', 
                           21: 'sequenceStatus', 22: 'VIndex', 23: 'n1Index',
                           24: 'n2Index', 25: 'DIndex', 26: 'JIndex'}
        self.productive = ("DataShare_1_Subject_Female_57yr_CD8_Memory\t"
                           + "UCSC-Kim-P01-01\tAATGCAACTTTAGCCACTCTGAAGATCCAG"
                           + "CCCTCAGAACCCAGGGACTCAGCTGTGTACTTCTGTGCCAGCGGTCT"
                           + "TCCACGACAGGGGGAAGAGACCCA\tCASGLPRQGEET\t"
                           + "3.5668468E-6\t15\t1.2044558E-6\t2\t39\t12\t"
                           + "(undefined)\tTRBV12-3, TRBV12-4\tTRBD1-1\t"
                           + "TRBJ2-5\t\t8\t2\t1\t3\t10\t0\tProductive\t63\t"
                           + "-1\t72\t82\t91")

    def test_parse_productive_clone(self):
        clone = adaptive.adaptive_parseline(self.productive, self.index2col)
        self.assertTrue(clone is not None)
        self.assertEqual(clone.count, 15)
        self.assertEqual(clone.freq, float("3.5668468E-6"))
        
        vgenes = ["TRBV12-3", "TRBV12-4"]
        self.assertTrue(set(clone.vgenes) == set(vgenes))
        jgenes = ["TRBJ2-5"]
        self.assertTrue(set(clone.jgenes) == set(jgenes))
        dgenes = ["TRBD1-1"]
        self.assertTrue(set(clone.dgenes) == set(dgenes))

        self.assertTrue(clone.productive)
        nuc = ("AATGCAACTTTAGCCACTCTGAAGATCCAGCCCTCAGAACCCAGGGACTCAGCTGTGTACT"
               + "TCTGTGCCAGCGGTCTTCCACGACAGGGGGAAGAGACC")
        self.assertEqual(clone.nuc, nuc)
        cdr3nuc = "TGTGCCAGCGGTCTTCCACGACAGGGGGAAGAGACC" 
        self.assertEqual(clone.cdr3nuc, cdr3nuc)
        aa = "NATLATLKIQPSEPRDSAVYFCASGLPRQGEET"
        self.assertEqual(clone.aa, aa)
        cdr3aa = "CASGLPRQGEET"
        self.assertEqual(clone.cdr3aa, cdr3aa)

        self.assertEqual(clone.lastvpos, 71)
        self.assertEqual(clone.firstdpos,  82)
        self.assertEqual(clone.lastdpos, 90)
        self.assertEqual(clone.firstjpos, 91)
         
        self.assertEqual(clone.vdel, 8)
        self.assertEqual(clone.d5del, 2)
        self.assertEqual(clone.d3del, 1)
        self.assertEqual(clone.jdel, 3)

class TestSequentaParseline(unittest.TestCase):
    '''Testing Sequenta parsing function
    '''
    def setUp(self):
        self.index2column = {0: 'Experiment_Name', 1: 'Lane', 2: 'Sample', 
                             3: 'Patient', 4: 'Clone_Index', 
                             5: 'Log10_Frequency', 6: 'Total_Read_Count',
                             7: 'Avg_QScore', 8: 'Read_Fraction_0_Errors',
                             9: 'J_Segment_Major_Gene', 
                             10: 'J_Segment_Major_Allele', 
                             11: 'J_Segment_Extension_Length', 
                             12: 'J_Segment_Deletion_Length', 
                             13: 'V_Segment_Major_Family', 
                             14: 'V_Segment_Major_Gene', 
                             15: 'V_Segment_Major_Allele', 
                             16: 'V_Segment_Extension_Length', 
                             17: 'V_Segment_Deletion_Length', 18: 'NDN_Length',
                             19: 'NDN_Effective_Length', 
                             20: 'N_Bases_adjacent_J',
                             21: 'N_Bases_adjacent_V', 
                             22: 'D_Segment_Major_Allele', 
                             23: 'D_Segment_length', 24: 'Is_Good_Frame', 
                             25: 'Clone_Sequence', 
                             26: 'Clone_V_Side_Extra_Sequence', 
                             27: 'CDR3_Sense_Sequence', 
                             28: 'Clone_Protein_Sequence'}
        self.productive = ("RNA_TCRB_AS_RInman_OldGroup_99\tRun036_6"
                           + "\tAS_RInman_OldGroup_99\tRinman_AS_99\t0"
                           + "\t-1.880\t16967\t35.98\t90.40\t"
                           + "TRBJ2-7\tTRBJ2-7*01\t45.00\t1.000\tTRBV7"
                           + "\tTRBV7-9\tTRBV7-9*01; TRBV7-9*02; TRBV7-9*03"
                           + "\t35.00\t3.000\t14.00\t12.00\t2\t2"
                           + "\tTRBD2*01; TRBD2*02\t10\ttrue\t"
                           + "TGTGACCGTGAGCCTGGTGCCCGGCCCGAAGTACTGCTCGTAGGAAG"
                           + "CGCTAGTCCCTGAGCTGCTGGCACAGAGATACATGGCCGAGTCCCCC"
                           + "\t\tTGTGCCAGCAGCTCAGGGACTAGCGCTTCCTACGAGC"
                           + "AGTACTTC\tGDSAMYLCASSSGTSASYEQYFGPGTRLTVT\n")

    def test_parse_productive_clone(self):
        clone = sequenta.sequenta_parseline(self.productive, self.index2column)
        self.assertTrue(clone is not None)
        self.assertEqual(clone.count, 16967)
        self.assertEqual(clone.freq, 10**(-1.880))
        
        vgenes = ["TRBV7-9"]
        self.assertTrue(set(clone.vgenes) == set(vgenes))
        valleles = ["TRBV7-9*01", "TRBV7-9*02", "TRBV7-9*03"]
        self.assertTrue(set(clone.valleles) == set(valleles))

        jgenes = ["TRBJ2-7"]
        self.assertTrue(set(clone.jgenes) == set(jgenes))
        jalleles = ["TRBJ2-7*01"]
        self.assertTrue(set(clone.jalleles) == set(jalleles))

        dalleles = ["TRBD2*01", "TRBD2*02"]
        self.assertTrue(set(clone.dalleles) == set(dalleles))
        dgenes = ["TRBD2"]
        self.assertTrue(set(clone.dgenes) == set(dgenes))

        self.assertTrue(clone.productive)
        nuc = ("GGGGACTCGGCCATGTATCTCTGTGCCAGCAGCTCAGGGACTAGCGCTTCCTACGAGCAGT"
               + "ACTTCGGGCCGGGCACCAGGCTCACGGTCACA")
        self.assertEqual(clone.nuc, nuc)
        cdr3nuc = "TGTGCCAGCAGCTCAGGGACTAGCGCTTCCTACGAGCAGTACTTC"
        self.assertEqual(clone.cdr3nuc, cdr3nuc)
        aa = "GDSAMYLCASSSGTSASYEQYFGPGTRLTVT"
        self.assertEqual(clone.aa, aa)
        cdr3aa = "CASSSGTSASYEQYF"
        self.assertEqual(clone.cdr3aa, cdr3aa)

        lastvpos = 35 - 2
        self.assertEqual(clone.lastvpos, lastvpos)

    def test_parse_water_clone(self):
        water = ("RNA_TCB_H2O_Micro_20130121_PCR289_Run047\tRun047_6"
                      + "H2O_Micro_20130121_PCR289\tWater\t5\t"
                      + "-1.754\t9\t32.95\t88.89\tTRBJ2-5\tTRBJ2-5*01"
                      + "43.00\t4.000\tTRBV15\tTRBV15\tTRBV15*02\t"
                      + "34.00\t.0\t17.00\t15.00\t9\t3\t"
                      + "TRBD2*02\t5\ttrue\tGAGCACCAGGAGCCGCGTGCC"
                      + "TGGCCCGAAGTACTGGGTCTCTCGCTCTCGACCCTCTCTGCTGGTGGCACAC"
                      + "AGGTACATGGCTGCGTCCCCC\t\tTGTGCCACCAGCAGAGAGGGT"
                      + "CGAGAGCGAGAGACCCAGTACTTC\tGDAAMYLCATSREGRERETQY"
                      + "FGPGTRLLVL\n")
        clone = sequenta.sequenta_parseline(water, self.index2column)
        self.assertTrue(clone is None)

    def test_parse_minfield_clone(self):
        clonestr = ("16967\t-1.880\tTGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC\t"
                    + "TRBV13\tTRBJ2-4")
        index2column = {1: 'Log10_Frequency', 0: 'Total_Read_Count',
                        2: 'Clone_Sequence', 3: 'V_Segment_Major_Gene',
                        4: 'J_Segment_Major_Gene'}
        clone = sequenta.sequenta_parseline(clonestr, index2column)
        self.assertTrue(clone is not None)

class TestProcessCloneInputs(unittest.TestCase):
    '''Testing see if the function reads relevant files. 
       Redudant names handling.
       Different input formats handling.
    '''
    def setUp(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        small_dir = "data/small"
        self.adaptive_dir = os.path.join(script_dir, small_dir, "adaptive")
        self.sequenta_dir = os.path.join(script_dir, small_dir, "sequenta")
        self.mitcr_dir = os.path.join(script_dir, small_dir, "mitcr")
        self.empty = os.path.join(script_dir, small_dir,
                                  "oddformats/emptyfile")
        self.wrongformat = os.path.join(script_dir, small_dir,
                                  "oddformats/wrongformatfile")
        self.minfield = os.path.join(script_dir, small_dir,
                                  "oddformats/minfield")
    
    def test_read_directory(self):
        # simplest case: 3 files, no extension
        name2sample = inputcommon.process_clone_inputs(self.sequenta_dir,
                                                       "sequenta")
        self.assertEqual(len(name2sample.keys()), 3)
        
        names = ["testsample1", "testsample2", "testsample3"]
        self.assertTrue(set(names) == set(name2sample.keys()))
        
        # wrong format type:
        self.assertRaises(inputcommon.UnrecognizedFormatError, 
                          inputcommon.process_clone_inputs,
                          self.sequenta_dir, "notsure")

        # test file extensions:
        name2sample = inputcommon.process_clone_inputs(self.adaptive_dir,
                                                       "adaptive", "tsv")
        self.assertEqual(len(name2sample.keys()), 2)
        
        names = ["female_57yr_cd8_memory_subset",
                 "male_35yr_cd8_memory_subset"]
        self.assertTrue(set(names) == set(name2sample.keys()))
        
        female_sample = name2sample["female_57yr_cd8_memory_subset"]
        numclones = len(female_sample.clones)
        self.assertEqual(numclones, 100)

        # test repetitive sample names:
        name2sample = inputcommon.process_clone_inputs(self.adaptive_dir,
                                                       "adaptive")
        self.assertEqual(len(name2sample.keys()), 3)
        
        names = ["female_57yr_cd8_memory_subset",
                 "male_35yr_cd8_memory_subset", "test"]
        self.assertTrue(set(names) == set(name2sample.keys()))
        
        female_sample = name2sample["female_57yr_cd8_memory_subset"]
        numclones = len(female_sample.clones)
        self.assertEqual(numclones, 109)

    def test_wrong_clone_file_format_handling(self):
        # test empty file
        self.assertRaises(inputcommon.FormatError, 
                          inputcommon.read_clone_file, self.empty)
        
        # test wrong format file
        self.assertRaises(inputcommon.FormatError, 
                          inputcommon.read_clone_file, self.wrongformat)
       
        # test minimum input file
        sample = inputcommon.read_clone_file(self.minfield)
        self.assertEqual(len(sample.clones), 10)

class TestReadGroupInfoFile(unittest.TestCase):
    '''Test reading the group info file
    '''
    def setUp(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        indir = "data/small/groupinfo"
        self.indir = os.path.join(script_dir, indir)

    def test_read_group_info_file(self):
        wrongfiles = ["empty.txt", "wrongmatchedline.txt", 
                      "zerogroupline.txt", "matchedwrong.txt"]
        for file in wrongfiles:
            filepath = os.path.join(self.indir, file)
            self.assertRaises(inputcommon.FormatError, 
                              inputcommon.read_group_info, filepath)
        
        nomatched = os.path.join(self.indir, "nomatched.txt")
        grs, group2samples, is_matched = inputcommon.read_group_info(nomatched)
        g2s = {'group1': ['s1.1', 's1.2', 's1.3'], 'group2': ['s2.1'],
               'group3': ['s3.1', 's3.2', 's3.3', 's3.4']}
        self.assertTrue(group2samples, g2s)
        self.assertTrue(not is_matched)
        
        matched = os.path.join(self.indir, "matchedcorrect.txt")
        grs, group2samples, is_matched = inputcommon.read_group_info(matched)
        g2s = {'group1': ['s1.1', 's1.2', 's1.3'], 
               'group2': ['s2.1', 's2.2', 's2.3'],
               'group3': ['s3.1', 's3.2', 's3.3']}
        self.assertTrue(group2samples, g2s)
        self.assertTrue(is_matched)
        
        rep = os.path.join(self.indir, "repetitivenames.txt")
        grs, group2samples, is_matched = inputcommon.read_group_info(matched)
        self.assertTrue(group2samples, g2s)
        self.assertTrue(is_matched)
        


if __name__ == '__main__':
    unittest.main()


