#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Parse MiTCR output file (mitcr.milaboratory.com)
'''

import os
import sys
import re

from aimseqtk.lib.clone import Clone
import aimseqtk.lib.common as libcommon


def mitcr_columns():
    cols = ["Read count", "Percentage", "CDR3 nucleotide sequence", 
            "CDR3 nucleotide quality", "Min quality", 
            "CDR3 amino acid sequence", "V alleles", "V segments", 
            "J alleles", "J segments", "D alleles", "D segments", 
            "Last V nucleotide position ", "First D nucleotide position",
            "Last D nucleotide position", "First J nucleotide position", 
            "VD insertions", "DJ insertions", "Total insertions"]
    return cols

def mitcr_parseline(line, index2col):
    items = line.strip().split('\t')
    if len(items) != len(index2col):
        sys.stderr.write("Incosistent number of columns between the following\
                          line and the header line, skipped it:\n\
                          Line:\n%s\n" %line)
        return None
    
    col2val = {}
    valid_cols = mitcr_columns()
    for i, col in enumerate(index2col):
        if col in valid_cols:
            col2val[col] = items[i]

    # Return None if line does not have minimum required fields.
    required_cols = ["Read count", "Percentage", "CDR3 nucleotide sequence",
                     "V segments", "J segments"]
    for c in required_cols:
        if c not in col2val or not col2val[c]:
            return None

    count = int(col2va['Read count'])
    freq = float(col2val['Percentage'])/100.0
    nuc = col2val['CDR3 nucleotide sequence']
    vgenes = col2val['V segments'].split(', ')
    jgenes = col2val['J segments'].split(', ')

    clone = Clone(count, freq, nuc, vgenes, jgenes, cdr3nuc=nuc)

    if 'D segments' in col2val:
        col2val.dgenes = col2val['D segments'].split(', ')
    if 'V alleles' in col2val:
        clone.valleles = col2val['V alleles'].split(', ')
    if 'J alleles' in col2val:
        clone.jalleles = col2val['J alleles'].split(', ')
    if 'D alleles' in col2val:
        clone.dalleles = col2val['D alleles'].split(', ')

    if 'CDR3 amino acid sequence' in col2val:
        clone.aa = col2val['CDR3 amino acid sequence']
    if 'Last V nucleotide position' in col2val:
        clone.lastvpos = int(col2val['Last V nucleotide position'])
    if 'First D nucleotide position' in col2val:
        clone.firstdpos = int(col2val['First D nucleotide position'])
    if 'Last D nucleotide position' in col2val:
        clone.lastdpos = int(col2val['Last D nucleotide position'])
    if 'First J nucleotide position' in col2val:
        clone.firstjpos = int(col2val['First J nucleotide position'])

    return clone

    

