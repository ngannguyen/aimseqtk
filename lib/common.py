#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys


def get_index2item(line):
    items = line.strip().split("\t")
    index2item = {}
    for i, item in enumerate(items):
        index2item[i] = item
    return index2item

def nt2aa(nt):
    # Translate nucleotide sequence to amino acid sequence
    codon2aa = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': 'Z', 'TGA': 'Z',
    'TTG': 'L', 'TCG': 'S', 'TAG': 'Z', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' 
    }
    aaseq = ''
    for i in xrange(0, len(nt)/3):
        codon = nt[i*3: i*3+3].upper()
        aaseq = aaseq + codon2aa[codon]
    return aaseq

def rc(nt):
    #Return a reverse complement of the input sequence.
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    rcnt = ''
    for i in xrange( len(nt) -1, -1, -1):
        if nt[i] in complement:
            rcnt += complement[nt[i]] 
        else:
            rcnt += nt[i]
    return rcnt


