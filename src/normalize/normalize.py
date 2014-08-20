#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Using Bioconductor's metagenomeSeq package normalization to normalize
the samples
'''

import os
import sys

import rpy2.robjects as robjs
import rpy2.robjects.numpy2ri as rpyn
from rpy2.rinterface import NARealType
import gzip
import cPickle as pickle

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.statcommon as statcommon
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon


def clone_matrix(colnames, clone2sample2size):
    # Convert into a matrix. Cols = Samples, Rows = Clones, Cells = Sizes
    # Sizes can be count or freq or normfreq (sizetype)
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
    return rows, rownames

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

def normalize_MRexp(rows, colnames, rownames):
    from rpy2.robjects.packages import importr
    mgs = importr("metagenomeSeq")
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

#========== job Objs ====
class CloneMatrixAgg(Target):
    def __init__(self, sams, indir, outfile):
        Target.__init__(self)
        self.sams = sams
        self.indir = indir
        self.outfile = outfile

    def run(self):
        objs = libcommon.load_pickledir(self.indir)
        colnames = self.sams
        rows = []
        rownames = []
        for obj in objs:
            rs, rns = clone_matrix(colnames, obj)
            rows.extend(rs)
            rownames.extend(rns)
        pickle.dump((rows, rownames, colnames), gzip.open(self.outfile, 'wb'))
        system("rm -Rf %s" % self.indir)

class CloneMatrix(Target):
    '''
    Convert into a matrix. Cols = Samples, Rows = Clones, Cells = Sizes
    Pickle matrix to outfile
    Sizes can be count or freq or normfreq (sizetype)
    '''
    def __init__(self, indir, outfile, workdir, sizetype):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
        self.workdir = workdir
        self.sizetype = sizetype

    def run(self):
        # get clone2sample2size
        sams = os.listdir(self.indir)
        self.addChildTarget(statcommon.GetClone2Samples(self.indir,
                                                self.workdir, self.sizetype))
        self.setFollowOnTarget(CloneMatrixAgg(sams, self.workdir,
                                                                 self.outfile))

class UpdateNormCount(Target):
    def __init__(self, infile, outfile, norm_col, names):
        Target.__init__(self)
        self.infile = infile
        self.outfile = outfile
        self.norm_col = norm_col
        self.names = names

    def run(self):
        name2nc = {}
        for i, name in enumerate(self.names):
            name2nc[name] = self.norm_col[i]
        
        clones = pickle.load(gzip.open(self.infile, "rb"))
        cloneids = self.names
        for clone in clones:
            normcount = 0.0
            id = clone.get_vseqj()
            assert id in cloneids
            #nc = self.norm_col.rx(id)
            nc = name2nc[id]
            if isinstance(nc, NARealType):
                nc = 0.0
            #clone.normcount = nc
            #HACK!!
            clone.count = nc
        pickle.dump(clones, gzip.open(self.outfile, "wb"))

class NormalizeMRexp2(Target):
    def __init__(self, matrix_file, samdir, outdir):
        Target.__init__(self)
        self.mx_file = matrix_file
        self.samdir = samdir
        self.outdir = outdir

    def run(self):
        self.logToMaster("NormalizeMRexp2")
        (rows, rownames, colnames) = pickle.load(gzip.open(self.mx_file, "rb"))
        #norm_matrix = normalize_MRexp(rows, rownames, colnames)
        norm_matrix = normalize_MRexp(rows, colnames, rownames)
        
        # get samples with normalized counts: 
        for sam in os.listdir(self.samdir):
            #norm_col = norm_matrix.rx[sam, True]
            norm_col = norm_matrix.rx[True, sam]
            #DEBUG
            #print sam
            #print len(norm_col)
            #print len(rownames)
            #rowindex = 0
            #for i, rn in enumerate(rownames):
            #    if rn == 'TCRBV12-03_CASSLGGVGAFF_TCRBJ01-01':
            #        rowindex = i
            #        print rowindex, rn
            #        break
            #print norm_col[rowindex]
            #col = []
            #index = -1
            #for i, colname in enumerate(colnames):
            #    if colname == sam:
            #        index = i
            #        print index, colname
            #        break
            ##print rows
            #print rows[rowindex*len(colnames) + index]
            #if sam == 'MH':            
            #    sys.exit(1)
            #END DEBUG
            samdir = os.path.join(self.samdir, sam)
            samout = os.path.join(self.outdir, sam)
            system("mkdir -p %s" % samout)
            for vj in os.listdir(samdir):
                vjin = os.path.join(samdir, vj)
                vjout = os.path.join(samout, vj)
                if vj == sam:  # sample file: out/sam/sam
                    system("cp %s %s" % (vjin, vjout))
                else:
                    assert len(rownames) == len(norm_col)
                    self.addChildTarget(UpdateNormCount(vjin, vjout, norm_col,
                                                        rownames))
        self.setFollowOnTarget(libcommon.CleanupFile(self.mx_file))

class NormalizeMRexp(Target):
    def __init__(self, indir, outdir, sizetype):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.sizetype = sizetype

    def run(self):
        matrix_tempdir = os.path.join(self.outdir, "matrix_temp") 
        matrix_file = os.path.join(self.outdir, "matrix_info.pickle")
        system("mkdir -p %s" % matrix_tempdir)
        self.addChildTarget(CloneMatrix(self.indir, matrix_file,
                                        matrix_tempdir, self.sizetype))
        self.setFollowOnTarget(NormalizeMRexp2(matrix_file, self.indir,
                                               self.outdir))


