#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Split clones by length
Indir:
    sample1/
        sample1, V1, V2, ...
    sample2/
        sample2, V1, V2, ...

Outdir:
    sample1/
        sample1
        V1/
            L1, L2, ...
        V2/
            L1, L2, ...
    sample2/
        ...
    ...
'''

import os
import sys
import numbers
import re
from math import log10, factorial
import cPickle as pickle
import gzip
import numpy as np
from scipy.stats import poisson

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon


def split_clones(clones):
    l2clones = {}
    for c in clones:
        l = len(c.aa)
        if l not in l2clones:
            l2clones[l] = [c]
        else:
            l2clones[l].append(c)
    return l2clones

class SplitClones(Target):
    def __init__(self, infile, outdir):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir

    def run(self):
        clones = pickle.load(gzip.open(self.infile, "rb"))
        l2clones = split_clones(clones)
        for l, lclones in l2clones.iteritems():
            outfile = os.path.join(self.outdir, str(l))
            pickle.dump(lclones, gzip.open(outfile, 'wb'))
        
class Setup(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            outsamdir = os.path.join(self.outdir, sam)
            system('mkdir -p %s' % outsamdir)
            for file in os.listdir(samdir):
                filepath = os.path.join(samdir, file)
                if file == sam:
                    system("cp %s %s" % (filepath, outsamdir))
                else:  # each V
                    outvdir = os.path.join(outsamdir, file)
                    system("mkdir -p %s" % outvdir)
                    self.addChildTarget(SplitClones(filepath, outvdir))

def main():
    usage = "%prog <indir> <outdir> [options]"
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    indir = args[0]
    outdir = args[1]
    i = Stack(Setup(indir, outdir)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.spl_clones_by_len import *
    main()

