#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Parse Adaptive Biotechnologies (http://www.adaptivebiotech.com/) data files
'''

import os
import sys

from aimseqtk.lib.clone import Clone
import aimseqtk.lib.common as libcommon


def adaptive_columns():
    cols = ['sequenceID', 'container', 'nucleotide', 'aminoAcid',
            'normalizedFrequency', 'normalizedCopy', 'rawFrequency',
            'copy', 'cdr3Length', 'VFamilyName', 'VGeneName', 'VTies',
            'DGeneName', 'JGeneName', 'JTies', 'VDeletion', 'd5Deletion',
            'd3Deletion', 'JDeletion', 'n2Insertion', 'n1Insertion',
            'sequenceStatus', 'VIndex', 'n1Index', 'n2Index', 'DIndex',
            'JIndex']
    return cols

def adaptive_parseline(line, index2col):
    items = line.strip().split('\t')
    if len(items) != len(index2col):
        sys.stderr.write("Incosistent number of columns between the following\
                          line and the header line, skipped it:\n\
                          Line:\n%s\n" %line)
        return None
    
    col2val = {}
    valid_cols = adaptive_columns()
    for i, col in enumerate(index2col):
        if col in valid_cols:
            col2val[col] = items[i]

    # Return None if line does not have minimum required fields.
    required_cols = ['normalizedCopy', 'normalizedFrequency', 'nucleotide',
                     'VGeneName', 'JGeneName']
    for c in required_cols:
        if c not in col2val or col2val[c] in ["(undefined)", ""]:
            return None

    count = int(col2val['normalizedCopy'])
    freq = float(col2val['normalizedFrequency'])
    nuc = col2val['nucleotide']
    vgenes = [col2val['VGeneName']]
    jgenes = [col2val['JGeneName']]

    # Clone with required fields
    clone = Clone(count, freq, nuc, vgenes, jgenes) 
    
    # Additional information if available
    # Gene info:
    if 'DGeneName' in col2val:
        if clone.dgenes[0] != '(undefined)':
            clone.dgenes = [col2val['DGeneName']]
    if 'VTies' in col2val:
        clone.vgenes = col2val['VTies'].split(", ")
    if 'JTies' in col2val:
        clone.jgenes = col2val['JTies'].split(", ")

    # Sequence ID, status and cdr3aa:
    if 'sequenceID' in col2val:
        clone.id = col2val['sequenceID']
    if 'sequenceStatus' in col2val:
        status = col2val['sequenceStatus']
        if status is not None and status.lower() == 'productive':
            clone.productive = True
        else:
            clone.productive = False
    if 'aminoAcid' in col2val:
        clone.cdr3aa = col2val['aminoAcid']

    # Junctional info:
    offset = 0
    if 'VIndex' in col2val:
        vindex = int(col2val['VIndex'])
        if clone.cdr3aa:
            cdr3len = clone.cdr3aa * 3
            endindex = max(len(clone.nuc), vindex + cdr3len)
            clone.cdr3nuc = clone.nuc[vindex: endindex]
        if clone.productive:
            # Make sure nuc is inframe:
            offset = vindex % 3
            nuclen = len(clone.nuc)
            endoffset = (nuclen - offset) % 3
            clone.nuc = clone.nuc[offset:  nuclen - endoffset]
            clone.aa = libcommon.nt2aa(clone.nuc)
    if 'DIndex' in col2val:
        clone.firstdpos = int(col2va['DIndex']) - offset
    if 'n2Index' in col2val:
        n2index = int(col2val['n2Index'])
        if n2index != -1:
            clone.lastvpos = n2index - 1 - offset
        elif clone.firstdpos:  # No d5ins
            clone.lastvpos = clone.firstdpos - 1 
    if 'JIndex' in col2val:
        clone.firstjpos = int(col2val['JIndex']) - offset
    if 'n1Index' in col2val:
        n1index = int(col2val['n1index'])
        if n1index != -1:
            clone.lastdpos = int(col2val['n1Index']) - 1 - offset
        elif clone.firstjpos:  # No d3ins
            clone.lastdpos = clone.firstjpos - 1

    # Deletion info:
    if 'VDeletion' in col2val:
        clone.vdel = int(col2val['VDeletion'])
    if 'd5Deletion' in col2val:
        clone.d5del = int(col2val['d5Deletion'])
    if 'd3Deletion' in col2val:
        clone.d3del = int(col2val['d3Deletion'])
    if 'JDeletion' in col2val:
        clone.jdel = int(col2val['JDeletion'])

    return clone

