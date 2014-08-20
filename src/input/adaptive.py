#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Parse Adaptive Biotechnologies (http://www.adaptivebiotech.com/) data files
Edit: Mon Jun 30 12:12:37 PDT 2014
Update to the current format of Adaptive TSV files (as of June 30 2014)
'''

import os
import sys

from aimseqtk.lib.clone import Clone
import aimseqtk.lib.common as libcommon

#nucleotide      aminoAcid       count   frequencyCount  cdr3Length      vMaxResolved    vFamilyName     vGeneName       vGeneAllele     vFamilyTies     vGeneNameTies   vGeneAlleleTies dMaxResolved    dFamilyName     dGeneName       dGeneAllele     dFamilyTies     dGeneNameTies   dGeneAlleleTies jMaxResolved    jFamilyName     jGeneName       jGeneAllele     jFamilyTies     jGeneNameTies   jGeneAlleleTies vDeletion       d5Deletion      d3Deletion      jDeletion       n2Insertion     n1Insertion     vIndex  n2Index dIndex  n1Index jIndex  estimatedNumberGenomes  sequenceStatus  cloneResolved   vOrphon dOrphon jOrphon vFunction       dFunction       jFunction       fractionNucleated       vAlignLength    vAlignSubstitutionCount vAlignSubstitutionIndexes       vAlignSubstitutionGeneThreePrimeIndexes

def adaptive_columns():
    #cols = ['nucleotide', 'aminoAcid',
    #        'count', 'frequencyCount'
    #        'cdr3Length', 'vFamilyName', 'vGeneName', 'VGeneNameTies',
    #        'dGeneName', 'jGeneName', 'jGeneNameTies', 'vDeletion', 'd5Deletion',
    #        'd3Deletion', 'jDeletion', 'n2Insertion', 'n1Insertion',
    #        'sequenceStatus', 'vIndex', 'n1Index', 'n2Index', 'dIndex',
    #        'jIndex']
    cols = ['nucleotide', 'aminoAcid', 'count', 'frequencyCount', 'cdr3Length',
            'vMaxResolved', 'vFamilyName', 'vGeneName', 'vGeneAllele',
            'vFamilyTies', 'vGeneNameTies', 'vGeneAlleleTies', 'dMaxResolved',
            'dFamilyName', 'dGeneName', 'dGeneAllele', 'dFamilyTies',
            'dGeneNameTies', 'dGeneAlleleTies', 'jMaxResolved', 'jFamilyName',
            'jGeneName', 'jGeneAllele', 'jFamilyTies', 'jGeneNameTies',
            'jGeneAlleleTies', 'vDeletion', 'd5Deletion', 'd3Deletion',
            'jDeletion', 'n2Insertion', 'n1Insertion', 'vIndex', 'n2Index',
            'dIndex', 'n1Index', 'jIndex', 'estimatedNumberGenomes',
            'sequenceStatus', 'cloneResolved', 'vOrphon', 'dOrphon', 'jOrphon',
            'vFunction', 'dFunction', 'jFunction', 'fractionNucleated',
            'vAlignLength', 'vAlignSubstitutionCount',
            'vAlignSubstitutionIndexes',
            'vAlignSubstitutionGeneThreePrimeIndexes']    
    return cols

def adaptive_parseline(line, index2col):
    items = line.strip('\n').split('\t')
    if len(items) != len(index2col):
        sys.stderr.write("Incosistent number of columns between the following\
                          line and the header line, skipped it:\n\
                          Line:\n%s\n" %line)
        return None
    
    col2val = {}
    valid_cols = adaptive_columns()
    for i, col in index2col.iteritems():
        if col in valid_cols:
            col2val[col] = items[i].replace("/", ", ")

    # Return None if line does not have minimum required fields.
    required_cols = ['count', 'frequencyCount', 'nucleotide',
                     'vGeneName', 'jGeneName', 'vGeneNameTies', 'jGeneNameTies']
    for c in required_cols:
        if c not in col2val:  # or col2val[c] in ["(undefined)", ""]:
            return None

    count = int(col2val['count'])
    freq = float(col2val['frequencyCount'])/100.0  # convert to non percentage
    nuc = col2val['nucleotide']
    vgene = col2val['vGeneName']
    if vgene == 'unresolved':
        vgenes = col2val['vGeneNameTies'].split(',')
    else:
        vgenes = [vgene]
    jgene = col2val['jGeneName']
    if jgene == 'unresolved':
        jgenes = col2val['jGeneNameTies'].split(',')
    else:
        jgenes = [jgene]


    # Clone with required fields
    clone = Clone(count, freq, nuc, vgenes, jgenes) 
    
    # Additional information if available
    # Gene info:
    if 'dGeneName' in col2val:
        dgenestr = col2val['dGeneName']
        if dgenestr == 'unresolved':
            clone.dgenes = col2val['dGeneNameTies'].split(',')
        else:
            clone.dgenes = [dgenestr]

    if 'sequenceStatus' in col2val:
        status = col2val['sequenceStatus'].lower()
        if status is not None and status == 'in':
            clone.productive = True
        elif status == 'out' or status == 'stop':
            clone.productive = False
        else:
            sys.stderr.write("Unknown status: %s\n" % status)
    if 'aminoAcid' in col2val:
        clone.cdr3aa = col2val['aminoAcid']

    # Junctional info:
    offset = 0
    if 'vIndex' in col2val:
        vindex = int(col2val['vIndex'])
        if clone.productive:
            # Make sure nuc is inframe:
            offset = vindex % 3
            nuclen = len(clone.nuc)
            endoffset = (nuclen - offset) % 3
            clone.nuc = clone.nuc[offset:  nuclen - endoffset]
            clone.aa = libcommon.nt2aa(clone.nuc)
        if clone.cdr3aa:
            cdr3len = len(clone.cdr3aa) * 3
            endindex = max(len(clone.nuc), vindex + cdr3len)
            clone.cdr3nuc = clone.nuc[vindex: endindex]
    if 'dIndex' in col2val:
        clone.firstdpos = int(col2val['dIndex']) - offset
    if 'n2Index' in col2val:
        n2index = int(col2val['n2Index'])
        if n2index != -1:
            clone.lastvpos = n2index - 1 - offset
        elif clone.firstdpos:  # No d5ins
            clone.lastvpos = clone.firstdpos - 1 
    if 'jIndex' in col2val:
        clone.firstjpos = int(col2val['jIndex']) - offset
    if 'n1Index' in col2val:
        n1index = int(col2val['n1Index'])
        if n1index != -1:
            clone.lastdpos = n1index - 1 - offset
        elif clone.firstjpos:  # No d3ins
            clone.lastdpos = clone.firstjpos - 1

    # Deletion info:
    if 'vDeletion' in col2val:
        clone.vdel = int(col2val['vDeletion'])
    if 'd5Deletion' in col2val:
        clone.d5del = int(col2val['d5Deletion'])
    if 'd3Deletion' in col2val:
        clone.d3del = int(col2val['d3Deletion'])
    if 'jDeletion' in col2val:
        clone.jdel = int(col2val['jDeletion'])

    return clone

