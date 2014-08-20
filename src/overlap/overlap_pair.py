#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Input: the db directory  (sample/V
       list of samples
       minimum number of samples
Output: a db directory that contains only clones that are shared by at
        least the specified number of samples
'''

import os
import sys
import cPickle as pickle
import gzip

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon
import aimseqtk.lib.statcommon as scommon
#import aimseqtk.src.overlap.overlap as overlap


def get_sample_obj_file(sam, outdir, db_dir):
    samdir = os.path.join(outdir, sam)
    system("mkdir -p %s" % samdir)
    samfile = os.path.join(db_dir, sam, sam)
    samobj = pickle.load(gzip.open(samfile, 'rb'))
    samobj.size = -1
    samobj.numclone = -1
    out_samfile = os.path.join(samdir, sam)
    pickle.dump(samobj, gzip.open(out_samfile, 'wb'))

def filter_clones(c2s, sams, minsam, maxsam):
    new_c2s = {}
    samset = set(sams)
    for c, sam2size in c2s.iteritems():
        clone_sams = set(sam2size.keys())
        common_sams = clone_sams.intersection(samset)
        if len(common_sams) >= minsam and (not maxsam or len(common_sams) <= maxsam):
            new_c2s[c] = list(common_sams)
    return new_c2s

class OverlapClones_Sample(Target):
    def __init__(self, infile, outfile, clones):
        Target.__init__(self)
        self.infile = infile
        self.outfile = outfile
        self.clones = clones

    def run(self):
        all_clones = pickle.load(gzip.open(self.infile, 'rb'))
        filter_clones = []
        for clone in all_clones:
            name = clone.get_vseqj()
            if name in self.clones:
                filter_clones.append(clone)
        pickle.dump(filter_clones, gzip.open(self.outfile, 'wb'))

class OverlapClones_V(Target):
    def __init__(self, db_dir, infile, v, minsam, maxsam, sams, outdir):
        Target.__init__(self)
        self.db_dir = db_dir
        self.infile = infile
        self.v = v
        self.minsam = minsam
        self.maxsam = maxsam
        self.sams = sams
        self.outdir = outdir

    def run(self):
        c2s = pickle.load(gzip.open(self.infile, 'rb'))
        filter_c2sams = filter_clones(c2s, self.sams, self.minsam, self.maxsam)
        sam2clones = lcommon.get_val2keys(filter_c2sams)
        for sam, clones in sam2clones.iteritems():
            infile = os.path.join(self.db_dir, sam, self.v)
            outfile = os.path.join(self.outdir, sam, self.v)
            self.addChildTarget(OverlapClones_Sample(infile, outfile, clones))

class ReadjustSample(Target):
    def __init__(self, indir):
        Target.__init__(self)
        self.indir = indir

    def run(self):
        sam = os.path.basename(self.indir)
        samfile = os.path.join(self.indir, sam)
        samobj = pickle.load(gzip.open(samfile, 'rb'))
        numclone = 0
        for v in os.listdir(self.indir):
            if v == sam:
                continue
            vfile = os.path.join(self.indir, v)
            clones = pickle.load(gzip.open(vfile, 'rb'))
            numclone += len(clones)
        samobj.numclone = numclone
        pickle.dump(samobj, gzip.open(samfile, 'wb'))

class ReadjustSamples(Target):
    def __init__(self, indir):
        Target.__init__(self)
        self.indir = indir

    def run(self):
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            self.addChildTarget(ReadjustSample(samdir))

class OverlapClones(Target):
    '''
    Record to outdir only clones that present in at least <minsam> number
    of input <sams>
    Outdir/
        sam1/
            sam1
            V1
            V2
            ...
        sam2/
        ...
    '''
    def __init__(self, db_dir, clone_dir, minsam, maxsam, sams, outdir):
        Target.__init__(self)
        self.db_dir = db_dir
        self.clone_dir = clone_dir
        self.minsam = minsam
        self.maxsam = maxsam
        self.sams = sams
        self.outdir = outdir

    def run(self):
        # Get the <sample> pickle file for each sample
        for sam in self.sams:
            get_sample_obj_file(sam, self.outdir, self.db_dir)
        for v in os.listdir(self.clone_dir):
            vfile = os.path.join(self.clone_dir, v)
            self.addChildTarget(OverlapClones_V(self.db_dir, vfile, v,
                    self.minsam, self.maxsam, self.sams, self.outdir))
        # Update the sample's <size> and <numclone> attrs
        self.setFollowOnTarget(ReadjustSamples(self.outdir))

class Setup(Target):
    def __init__(self, args, options):
        Target.__init__(self)
        self.args = args
        self.c2s_dir = options.c2s2s_dir
        self.maxsam = options.maxsam

    def run(self):
        db_dir = self.args[0]
        sam_file = self.args[1]
        minsam = int(self.args[2])
        outdir = self.args[3]

        if not os.path.exists(outdir):
            system("mkdir -p %s" % outdir)
        samples = lcommon.read_list(sam_file)
        
        c2s_dir = self.c2s_dir
        if c2s_dir is None:
            global_dir = self.getGlobalTempDir()
            c2s_dir = os.path.join(global_dir, "clone2sam2size_dir") 
            system("mkdir -p %s" % c2s_dir)
            sizetype = 'freq'
            self.addChildTarget(scommon.GetClone2Samples(db_dir, c2s_dir,
                                                         sizetype, samples))
        self.setFollowOnTarget(OverlapClones(db_dir, c2s_dir, minsam,
                                             self.maxsam, samples, outdir))

def main():
    usage = ("%prog <db_dir> <sample_list> <minimum_sample> <outdir>")
    parser = lcommon.init_options(usage)
    parser.add_option('-c', '--c2s2s_dir', default=None,
                      help=('Directory contains clone2sample2size info. ' +
                            'Default=%default.'))
    parser.add_option('--maxsam', dest='maxsam', type='int', default=None)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    i = Stack(Setup(args, options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.overlap.overlap_pair import *
    main()

