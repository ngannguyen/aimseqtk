#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Input: the db directory  (sample/V
Output: readjust the sample objs 
'''

import os
import sys
import cPickle as pickle
import gzip

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon

class AdjustSample(Target):
    def __init__(self, indir, sam):
        Target.__init__(self)
        self.indir = indir
        self.sam = sam

    def run(self):
        # Get the <sample> pickle file for each sample
        indir = os.path.join(self.indir, self.sam)
        samfile = os.path.join(indir, self.sam)
        samobj = pickle.load(gzip.open(samfile, 'rb'))
        numclone = 0
        for v in os.listdir(indir):
            if v == self.sam:
                continue
            vfile = os.path.join(indir, v)
            clones = pickle.load(gzip.open(vfile, 'rb'))
            numclone += len(clones)
        samobj.numclone = numclone
        pickle.dump(samobj, gzip.open(samfile, 'wb'))

class Setup(Target):
    def __init__(self, args):
        Target.__init__(self)
        self.indir = args[0]

    def run(self):
        for sam in os.listdir(self.indir):
            self.addChildTarget(AdjustSample(self.indir, sam))

def main():
    usage = ("%prog <db_dir>")
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    i = Stack(Setup(args)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.readjust_samples import *
    main()

