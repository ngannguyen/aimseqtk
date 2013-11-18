#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Parse Sequenta data files (sequentainc.com)
'''

import os
import sys
import re

from aimseqtk.lib.clone import Clone
import aimseqtk.lib.common as libcommon


def sequenta_columns():
    cols = ['Experiment_Name', 'Lane', 'Sample', 'Patient', 'Clone_Index', 
            'Log10_Frequency', 'Total_Read_Count', 'Avg_QScore', 
            'Read_Fraction_0_Errors', 'J_Segment_Major_Gene', 
            'J_Segment_Major_Allele', 'J_Segment_Extension_Length', 
            'J_Segment_Deletion_Length', 'V_Segment_Major_Family', 
            'V_Segment_Major_Gene', 'V_Segment_Major_Allele', 
            'V_Segment_Extension_Length', 'V_Segment_Deletion_Length', 
            'NDN_Length', 'NDN_Effective_Length', 'N_Bases_adjacent_J', 
            'N_Bases_adjacent_V', 'D_Segment_Major_Allele', 'D_Segment_length',
            'Is_Good_Frame', 'Clone_Sequence', 'Clone_V_Side_Extra_Sequence',
            'CDR3_Sense_Sequence', 'Clone_Protein_Sequence']
    return cols

def sequenta_getaa(nuc):
    if re.match("N", nuc):
        assert len(nuc) >= 3
        aa = "X" + libcommon.nt2aa(nuc[3:])
    else:
        aa = libcommon.nt2aa[nuc]

def sequenta_parseline(line, index2col):
    items = line.strip().split('\t')
    if len(items) != len(index2col):
        sys.stderr.write("Incosistent number of columns between the following\
                          line and the header line, skipped it:\n\
                          Line:\n%s\n" %line)
        return None
    
    col2val = {}
    valid_cols = sequenta_columns()
    for i, col in enumerate(index2col):
        if col in valid_cols:
            col2val[col] = items[i]

    # Return None if line does not have minimum required fields.
    required_cols = ['Total_Read_Count', 'Log10_Frequency', 'Clone_Sequence',
                     'V_Segment_Major_Gene', 'J_Segment_Major_Gene'] 
    for c in required_cols:
        if c not in col2val or col2val[c] in ['NAN', '', '-']:
            return None

    count = int(col2va['Total_Read_Count'])
    freq = 10**float(col2val['Log10_Frequency'])
    nuc = col2val['Clone_Sequence']
    vgenes = col2val['V_Segment_Major_Gene'].split('; ')
    jgenes = col2val['J_Segment_Major_Gene'].split('; ')
    # Clone with required fields
    clone = Clone(count, freq, nuc, vgenes, jgenes) 

    # Additional information if available
    # Gene info:
    dstr = col2val['D_Segment_Major_Allele']
    if dstr not in ['NAN', '', '-']:
        dalleles = dstr.split('; ')
        dgenes = []
        for d in dalleles:
            dgene = d.split("*")[0]
            if dgene not in dgenes:
                dgenes.append(dgene)
        clone.dgenes = dgenes
        clone.dalleles = dalleles
    if 'V_Segment_Major_Allele' in col2val:
        clone.valleles = col2val['V_Segment_Major_Allele'].split('; ')
    if 'J_Segment_Major_Allele' in col2val:
        clone.jalleles = col2val['J_Segment_Major_Allele'].split('; ')
      
    # Sequence ID, status and cdr3aa:
    if 'Sample' in col2val:
        clone.samplename = col2val['Sample']
    if 'Patient' in col2val:
        clone.patient = col2val['Patient']
    if 'Clone_Index' in col2val:
        clone.id = col2val['Clone_Index']
    if 'Is_Good_Frame' in col2val:
        if col2val['Is_Good_Frame'].lower() == 'true':
            clone.productive = True
        else:
            clone.productive = False
    if 'Clone_Protein_Sequence' in col2val:
        clone.aa = col2val['Clone_Protein_Sequence'].replace('*', 'Z')
    
    if 'CDR3_Sense_Sequence' in col2val:
        clone.cdr3nuc = col2val['CDR3_Sense_Sequence']
    if not re.search(clone.cdr3nuc, clone.nuc):
        clone.nuc = libcommom.rc(clone.nuc)
    try:
        sequenta_getaa(clone.cdr3nuc)
    except:  # return None if cannot translate cdr3nuc
        return None

    # Make sure nuc is in frame
    cdr3start = re.search(clone.cdr3nuc, clone.nuc).start()
    offset = cdr3start % 3
    nuclen = len(clone.nuc)
    endoffset = (nuclen - offset) % 3
    clone.nuc = clone.nuc[offset: nuclen - endoffset]
   
    # Junctional info:
    if 'V_Segment_Extension_Length' in col2val:
        vins = int(col2val['V_Segment_Extension_Length'])
        clone.lastvpos = vins - 1 - offset
        if 'N_Bases_adjacent_V' in col2val:
            d5ins = col2val['N_Bases_adjacent_V']
            if not d5ins.startswith('-') and d5ins not in ['', 'NAN']:
                clone.firstdpos = clone.lastvpost + int(d5ins)
    if 'J_Segment_Extension_Length' in col2val:
        jins = int(col2val['J_Segment_Extension_Length']
        clone.firstjpos = len(clone.nuc) - jins
        if 'N_Bases_adjacent_J' in col2val:
            d3ins = col2val['N_Bases_adjacent_J']
            if not d3ins.startswith('-') and d3ins not in ['', 'NAN']:
                clone.lastdpos = clone.firstjpost - int(d3ins) - 1

    # Deletions:
    if 'V_Segment_Deletion_Length' in col2val:
        vdel = col2val['V_Segment_Deletion_Length']
        if not vdel.startswith('-') and vdel not in ['', 'NAN']:
            clone.vdel = int(vdel)
    if 'J_Segment_Deletion_Length' in col2val:
        jdel = col2val['J_Segment_Deletion_Length']
        if not jdel.startswith('-') and jdel not in ['', 'NAN']:
            clone.jdel = int(jdel)
    
    return clone
