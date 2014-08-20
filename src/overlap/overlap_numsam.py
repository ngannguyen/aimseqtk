#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
1/ split the clones by number of samples that share them
2/ randomly select n clones for each number of samples
3/ print summary results: numsam2clonecount

Input: the clone2sample2size directory
       minimum number of samples
Output:
    outdir/
        clones/
            2 (file of all clones shares by 2 samples)
            3 (file of all clones shares by 3 samples)
            ...
        sampling_clones.txt
        numsam2clonecount.txt
'''


import os
import sys
import random
import cPickle as pickle
import gzip

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon
import aimseqtk.lib.statcommon as scommon


def split_clones_by_numsam(clone2sams, minsam, maxsam):
    numsam2clones = {}
    for c, sams in clone2sams.iteritems():
        numsam = len(sams)
        if numsam >= minsam:
            if maxsam is not None and numsam > maxsam:
                continue
            if numsam not in numsam2clones:
                numsam2clones[numsam] = [c]
            else:
                numsam2clones[numsam].append(c)
    return numsam2clones
    
class SplitClonesByNumsam(Target):
    def __init__(self, infile, outdir, minsam, maxsam):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.minsam = minsam
        self.maxsam = maxsam

    def run(self):
        clone2sams = pickle.load(gzip.open(self.infile, 'rb'))
        numsam2clones = split_clones_by_numsam(clone2sams, self.minsam,
                                               self.maxsam)
        for num, clones in numsam2clones.iteritems():
            outfile = os.path.join(self.outdir, str(num))
            pickle.dump(clones, gzip.open(outfile, "wb"))

class SplitClonesByNumsamAgg(Target):
    def __init__(self, indir, outdir, sampling):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.sampling = sampling

    def run(self):
        nums = []
        for v in os.listdir(self.indir):
            vdir = os.path.join(self.indir, v)
            for num in os.listdir(vdir):
                if num not in nums:
                    nums.append(num)

        numsam2clonecount = {}
        f2 = None
        if self.sampling is not None:
            samplingfile = os.path.join(self.outdir, "subclones.txt")
            f2 = open(samplingfile, 'w')
            f2.write("#Clone\tNum_samples\n")

        clonedir = os.path.join(self.outdir, "clones")
        system("mkdir -p %s" % clonedir)

        for num in nums:
            outfile = os.path.join(clonedir, num)
            clones = []
            for v in os.listdir(self.indir):
                vfile = os.path.join(self.indir, v, num)
                if os.path.exists(vfile):
                    vclones = pickle.load(gzip.open(vfile, 'rb'))
                    clones.extend(vclones)
            f = open(outfile, 'w')
            for c in clones:
                f.write("%s\n" % c)
            f.close()
            
            numsam2clonecount[num] = len(clones)
            if self.sampling is not None:
                if self.sampling > len(clones):
                    #subclones = clones
                    subclones = []
                else:
                    subclones = random.sample(clones, self.sampling)
                for sc in subclones:
                    f2.write("%s\t%s\n" % (sc, num))
        f.close()
        if f2:
            f2.close()
        sumfile = os.path.join(self.outdir, "numsam2clonecount.txt")
        f3 = open(sumfile, 'w')
        f3.write("#Numsam\tClonecount\n")
        nums = sorted([int(n) for n in nums])
        for n in nums:
            n = str(n)
            clonecount = numsam2clonecount[n]
            f3.write("%s\t%d\n" % (n, clonecount))
        f3.close()

class SplitClones(Target):
    def __init__(self, indir, minsam, maxsam, outdir, sampling):
        Target.__init__(self)
        self.indir = indir
        self.minsam = minsam
        self.maxsam = maxsam
        self.outdir = outdir
        self.sampling = sampling

    def run(self):
        global_dir = self.getGlobalTempDir()
        for v in os.listdir(self.indir):
            vfile = os.path.join(self.indir, v)
            vdir = os.path.join(global_dir, v)
            system('mkdir -p %s' % vdir)
            self.addChildTarget(SplitClonesByNumsam(vfile, vdir,
                                                    self.minsam, self.maxsam))
        self.setFollowOnTarget(SplitClonesByNumsamAgg(global_dir, self.outdir,
                                                      self.sampling))

class Setup(Target):
    def __init__(self, args, options):
        Target.__init__(self)
        self.indir = args[0]
        self.minsam = int(args[1])
        self.outdir = args[2]
        self.sampling = options.sampling
        self.maxsam = options.maxsam
        self.db = options.db

    def run(self):
        c2sdir = self.indir
        if self.db:
            c2sdir = os.path.join(self.outdir, "clone2sam2size_db")
            system("mkdir -p %s" % c2sdir)
            sizetype = 'freq'
            self.addChildTarget(scommon.GetClone2Samples(self.indir, c2sdir,
                                                         sizetype))
        self.setFollowOnTarget(SplitClones(c2sdir, self.minsam, self.maxsam,
                                           self.outdir, self.sampling))

def main():
    usage = ("%prog <indir> <minimum_sample> <outdir>")
    parser = lcommon.init_options(usage)
    parser.add_option('-m', '--maxsam', dest='maxsam', type='int',
                      default=None, help=('max number of samples'))
    parser.add_option('--db', dest='db', action='store_true', default=False,
                      help=('Specified if indir is the db dir instead of the' +
                            ' clone2sample2size dir. Default=%default'))
    parser.add_option('-s', '--sampling', dest='sampling', default=None,
                      type='int', help=('number of clones in each numshare ' +
                                        'category to take. default=%default'))
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    i = Stack(Setup(args, options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.overlap.overlap_numsam import *
    main()

