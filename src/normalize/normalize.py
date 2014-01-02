#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Using Bioconductor's metagenomeSeq package normalization to normalize
the samples
'''

import os
import sys
import gzip
import cPickle as pickle
from optparse import OptionGroup

import rpy2.robjects as robjs
import rpy2.robjects.numpy2ri as rpyn
from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
import aimseqtk.lib.statcommon as statcommon


def clone_matrix(samples, sizetype):
    # Convert into a matrix. Cols = Samples, Rows = Clones, Cells = Sizes
    # Sizes can be count or freq or normfreq (sizetype)
    clone2sample2size = statcommon.get_clone2samples(samples, sizetype)
    colnames = [s.name for s in samples]  # sample names
    rownames = clone2sample2size.keys()  # clone ids: v_cdr3_j
    rows = []
    for clone, sam2size in clone2sample2size.iteritems():
        row = []
        for sam in colnames:
            size = 0
            if sam in sam2size:
                size = sam2size[sam]
            row.append(size)
        rows.extend(row)
    return rows, colnames, rownames

def get_R_matrix(rows, colnames, rownames):
    numcol = len(colnames)
    numrow = len(rownames)
    v = robjs.FloatVector(rows)
    m = robjs.r['matrix'](v, ncol=numcol, nrow=numrow, byrow=True,
                          dimnames=[rownames, colnames])
    return m

def get_meta_matrix(group2samples):
    # Metadata information: rows = samples; cols = group
    colnames = ['Group']
    rownames = []
    rows = []
    for group, samples in group2samples.iteritems():
        rownames.extend(samples)
        rows.extend([group] * len(samples))
    return rows, colnames, rownames

#def prepare_MRexp(samples, sizetype, group2samples):
    ## get phenotype matrix:
    #r, cn, rn = get_meta_matrix(group2samples)
    #pheno_matrix = get_R_matrix(r, cn, rn)

def normalize_MRexp(samples, sizetype):
    from rpy2.robjects.packages import importr
    mgs = importr("metagenomeSeq")
    # get count matrix:
    rows, colnames, rownames = clone_matrix(samples, sizetype)
    count_matrix = get_R_matrix(rows, colnames, rownames)
    # prepare MRexperiment object:
    mrexp = mgs.newMRexperiment(count_matrix)
    # normalized using CSS:
    #normstat = mgs.cumNormStat(mrexp)
    #mrexp2 = mgs.cumNorm(mrexp, p=normstat)
    norm_count_matrix = mgs.MRcounts(mrexp, norm=True)
    return norm_count_matrix

def matrix_to_normcount(matrix, samples):
    # read normalized count from matrix (row=clone; col=sample) and
    # update the samples' clones
    samplenames = matrix.colnames
    cloneids = matrix.rownames
    for sample in samples:
        if sample.name not in samplenames:
            sys.stderr.write(("Warning: sample %s does not have normalized "
                              % sample.name + "count."))
            continue
        for clone in sample.clones:
            normcount = 0.0
            for id in clone.vjseq_ids:
                assert id in cloneids
                nc = matrix.rx[id, sample.name]
                normcount = normcount + rpyn.ri2numpy(nc)[0]
            clone.set_normcount(normcount)
    return samples


