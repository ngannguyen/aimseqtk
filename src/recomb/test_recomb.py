'''Compute the likelihood of generating sequences of each specific length
'''

import os
import sys
from math import log10
import cPickle as pickle
import gzip

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon
import aimseqtk.src.recomb.recomb_common as rcommon


class GetLlh_VJ(Target):
    def __init__(self, model, infile, outfile):
        Target.__init__(self)
        self.model = model
        self.infile = infile
        self.outfile = outfile

    def run(self):
        llhs = []
        jclones = pickle.load(gzip.open(self.infile, 'rb'))
        for c in jclones:
            clonellh = rcommon.ntclone_likelihood(c, self.model)
            llhs.append(clonellh)
        sumllh = sum([10 ** llh for llh in llhs])
        pickle.dump(sumllh, gzip.open(self.outfile, 'wb')) 

class GetLlh_V_Cleanup(Target):
    def __init__(self, basename):
        Target.__init__(self)
        self.basename = basename

    def run(self):
        dbdir = "%s-db_jsplit" % self.basename
        system("rm -Rf %s" % dbdir)
        llhdir = "%s-llh_jsplit" % self.basename
        system("rm -Rf %s" % llhdir)

class GetLlh_VJ_Agg(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        sumllh = 0.0
        for j in os.listdir(self.indir):
            jfile = os.path.join(self.indir, j)
            llh = pickle.load(gzip.open(jfile, 'rb'))
            sumllh += llh
        pickle.dump(sumllh, gzip.open(self.outfile, 'wb'))
        self.setFollowOnTarget(GetLlh_V_Cleanup(self.outfile))

class GetLlhs(Target):
    def __init__(self, db_file, model, outfile):
        Target.__init__(self)
        self.model = model
        self.db_file = db_file
        self.outfile = outfile
    
    def run(self):
        # split db clones by j
        js = []
        tempdir = "%s-db_jsplit" % self.outfile
        system("mkdir -p %s" % tempdir)
        clones = pickle.load(gzip.open(self.db_file, "rb"))
        j2clones = lcommon.split_clones_by_j(clones)
        for j, jclones in j2clones.iteritems():
            tempjfile = os.path.join(tempdir, j)  #  v/l-db_jsplit/jfile
            pickle.dump(jclones, gzip.open(tempjfile, 'wb'))
            if j not in js:
                js.append(j)
        self.logToMaster("Done spliting clones by j for %s\n" % self.outfile)
        
        tempoutdir = "%s-llh_jsplit" % self.outfile  #  v/l-llh_jsplit/jllh
        system("mkdir -p %s" % tempoutdir)
        for j in js:
            jinfile = os.path.join(tempdir, j)
            joutfile = os.path.join(tempoutdir, j)
            self.addChildTarget(GetLlh_VJ(self.model, jinfile, joutfile))
        self.setFollowOnTarget(GetLlh_VJ_Agg(tempoutdir, self.outfile)) 

class GetLlhsAgg(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        llh = 0.0
        for v in os.listdir(self.indir):
            vdir = os.path.join(self.indir, v)
            for l in os.listdir(vdir):
                lfile = os.path.join(vdir, l)
                vl_llh = pickle.load(gzip.open(lfile, "rb"))
                llh += vl_llh
        f = open(self.outfile, 'w')
        f.write("%f\n" % llh)
        f.close()

class Setup(Target):
    def __init__(self, db_dir, model_file, outfile, options):
        Target.__init__(self)
        self.db_dir = db_dir
        self.model_file = model_file
        self.outfile = outfile
        self.options = options

    def run(self):
        model = pickle.load(gzip.open(self.model_file, 'rb'))
        self.logToMaster("Done loading model.\n")
        
        global_dir = self.getGlobalTempDir()
        for v in os.listdir(self.db_dir):
            if v == os.path.basename(self.db_dir):
                continue
            vdir = os.path.join(self.db_dir, v)
            vdirout = os.path.join(global_dir, v)
            system('mkdir -p %s' % vdirout)
            for l in os.listdir(vdir):
                lfile = os.path.join(vdir, l)
                lfileout = os.path.join(vdirout, l)  # global_dir/v/l
                self.addChildTarget(GetLlhs(lfile, model, lfileout))
        self.setFollowOnTarget(GetLlhsAgg(global_dir, self.outfile))

def main():
    usage = "%prog <db_dir> <model_file> <outfile> [options]"
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    db_dir = args[0]
    model_file = args[1]
    outfile = args[2]

    i = Stack(Setup(db_dir, model_file, outfile, options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.test_recomb import *
    main()

