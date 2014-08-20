#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Parse Sequenta data files (sequentainc.com)
'''

import os
import sys
import re
import random

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
        aa = libcommon.nt2aa(nuc)
    return aa

def get_j_groups(jgenes):
    # return the groups (1 or 2 or both) that the jgenes are in
    jgroups = []
    for jgene in jgenes:
        group = jgene.lstrip('TRBJ').split('-')[0]
        if group not in jgroups:
            jgroups.append(group)
    return jgroups

def get_ddels(ddel):
    # randomly select a combination of (d5del, d3del) out of all
    # possible combinations given ddel
    if ddel <= 0:
        return 0, 0
    else:
        try:
            d5del = random.randint(0, ddel) 
            d3del = ddel - d5del
            return d5del, d3del
        except:
            print ddel
            sys.exit()

def sequenta_parseline(line, index2col):
    items = line.strip("\n").split('\t')
    if len(items) != len(index2col):
        sys.stderr.write("Incosistent number of columns between the following\
                          line and the header line, skipped it:\n\
                          Line:\n%s\n" %line)
        return None
   
    col2val = {}
    valid_cols = sequenta_columns()
    for i, col in index2col.iteritems():
        if col in valid_cols:
            col2val[col] = items[i]

    # Return None if clone is "Water"
    if 'Patient' in col2val and col2val['Patient'] == 'Water':
        return None

    # Return None if line does not have minimum required fields.
    required_cols = ['Total_Read_Count', 'Log10_Frequency', 'Clone_Sequence',
                     'V_Segment_Major_Gene', 'J_Segment_Major_Gene'] 
    for c in required_cols:
        if c not in col2val or col2val[c] in ['NAN', '', '-']:
            return None

    count = libcommon.soft_int(col2val['Total_Read_Count'])
    try:
        freq = 10**float(col2val['Log10_Frequency'])
    except:  # Return None if clone does not have a valid frequency
        return None
    nuc = col2val['Clone_Sequence']
    vgenes = col2val['V_Segment_Major_Gene'].split('; ')
    jgenes = col2val['J_Segment_Major_Gene'].split('; ')
    # Clone with required fields
    clone = Clone(count, freq, nuc, vgenes, jgenes) 
    
    # Additional information if available
    # Gene info:
    if 'D_Segment_Major_Allele' in col2val:
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
    if not clone.dgenes:  # no dgenes info
        jgroups = get_j_groups(clone.jgenes)
        if ['1'] == jgroups:
            clone.dgenes = ['TRBD1']
        else:
            clone.dgenes = [random.choice(['TRBD1', 'TRBD2'])]

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
    
    offset = 0
    if 'CDR3_Sense_Sequence' in col2val:
        clone.cdr3nuc = col2val['CDR3_Sense_Sequence']
        if not re.search(clone.cdr3nuc, clone.nuc):
            clone.nuc = libcommon.rc(clone.nuc)
        try:
            cdr3aa = sequenta_getaa(clone.cdr3nuc)
            clone.cdr3aa = cdr3aa
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
        vins = libcommon.soft_int(col2val['V_Segment_Extension_Length'])
        clone.lastvpos = vins - 1 - offset
        if 'N_Bases_adjacent_V' in col2val:
            d5ins = col2val['N_Bases_adjacent_V']
            if not d5ins.startswith('-') and d5ins not in ['', 'NAN']:
                clone.firstdpos = clone.lastvpos + int(d5ins) + 1
    if 'J_Segment_Extension_Length' in col2val:
        jins = libcommon.soft_int(col2val['J_Segment_Extension_Length'])
        clone.firstjpos = len(clone.nuc) - jins
        if 'N_Bases_adjacent_J' in col2val:
            d3ins = col2val['N_Bases_adjacent_J']
            if not d3ins.startswith('-') and d3ins not in ['', 'NAN']:
                clone.lastdpos = clone.firstjpos - int(d3ins) - 1

    # Deletions:
    if 'V_Segment_Deletion_Length' in col2val:
        vdel = col2val['V_Segment_Deletion_Length']
        if not vdel.startswith('-') and vdel not in ['', 'NAN']:
            clone.vdel = libcommon.soft_int(vdel)
    if 'J_Segment_Deletion_Length' in col2val:
        jdel = col2val['J_Segment_Deletion_Length']
        if not jdel.startswith('-') and jdel not in ['', 'NAN']:
            clone.jdel = libcommon.soft_int(jdel)
    
    # Special treatment for D info:
    d2fulllen = {'TRBD1': 12, 'TRBD2': 16}
    if 'D_Segment_length' in col2val:
        dgene = clone.dgenes[0]
        dfulllen = d2fulllen[dgene]
        dlen = col2val['D_Segment_length']
        if not dlen.startswith('-') and dlen not in ['', 'NAN']:
            ddel = dfulllen - int(dlen)
            clone.d5del, clone.d3del = get_ddels(ddel)
            #clone.d5del = ddel / 2.0
            #clone.d3del = ddel - clone.d5del
        else:  # all D was deleted
            clone.d5del, clone.d3del = get_ddels(dfulllen)
            #clone.d5del = dfulllen / 2
            #clone.d3del = dfulllen - clone.d5del
            ndn = clone.firstjpos - clone.lastvpos
            clone.firstdpos = clone.lastvpos + ndn/2 + 1
            clone.lastdpos = clone.firstdpos - 1

    return clone
