
import os
import sys
import gzip
import cPickle as pk

from sonLib.bioio import system 
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

import aimseqtk.lib.common as lcommon


class GetLenClones(Target):
    def __init__(self, indir, outfile, minlen, maxlen):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
        self.minlen = minlen
        self.maxlen = maxlen

    def run(self):
        vclones = []
        for l in os.listdir(self.indir):
            lfile = os.path.join(self.indir, l)
            l = int(l)
            if self.minlen <= l and l <= self.maxlen:
                clones = pk.load(gzip.open(lfile, 'rb'))
                vclones.extend(clones)
        pk.dump(vclones, gzip.open(self.outfile, 'wb'))

class Setup(Target):
    def __init__(self, args):
        Target.__init__(self)
        self.indir = args[0]
        self.minlen = int(args[1])
        self.maxlen = int(args[2])
        self.outdir = args[3]

    def run(self):
        system('mkdir -p %s' % self.outdir)
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            outsamdir = os.path.join(self.outdir, sam)
            system('mkdir -p %s' % outsamdir)
            for v in os.listdir(samdir):
                vpath = os.path.abspath(os.path.join(samdir, v))
                outfile = os.path.abspath(os.path.join(outsamdir, v))
                if v == sam:
                    system("ln -s %s %s" % (vpath, outfile))
                else:
                    self.addChildTarget(GetLenClones(vpath, outfile,
                                           self.minlen, self.maxlen))

def main():
    parser = lcommon.init_options() 
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    i = Stack(Setup(args)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.get_db_len import *
    main()
