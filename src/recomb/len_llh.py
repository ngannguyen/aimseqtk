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


def get_union_children(indir):
    children = []
    for s in os.listdir(indir):
        for file in os.listdir(os.path.join(indir, s)):
            if file == s:
                continue
            if file not in children:
                children.append(file)
    return children

def get_lens(db_dir):
    lens = []
    for s in os.listdir(db_dir):
        sdir = os.path.join(db_dir, s)
        for v in os.listdir(sdir):
            if v == s:
                continue
            vdir = os.path.join(sdir, v)
            for l in os.listdir(vdir):
                if l not in lens:
                    lens.append(l)
    return lens

class GetLenLlh_VJ(Target):
    def __init__(self, model, indir, j, outfile):
        Target.__init__(self)
        self.model = model
        self.indir = indir
        self.j = j
        self.outfile = outfile

    def run(self):
        llhs = []
        events = []
        for sam in os.listdir(self.indir):
            jfile = os.path.join(self.indir, sam, self.j)
            if os.path.exists(jfile):
                jclones = pickle.load(gzip.open(jfile, 'rb'))
                for c in jclones:
                    if not rcommon.visited_event(events, c):
                        events.append(c)
                        clonellh = rcommon.ntclone_likelihood(c, self.model)
                        llhs.append(clonellh)
        sumllh = sum([10 ** llh for llh in llhs])
        pickle.dump(sumllh, gzip.open(self.outfile, 'wb')) 

class GetLenLlh_V_Cleanup(Target):
    def __init__(self, basename):
        Target.__init__(self)
        self.basename = basename

    def run(self):
        dbdir = "%s-db_jsplit" % self.basename
        system("rm -Rf %s" % dbdir)
        llhdir = "%s-llh_jsplit" % self.basename
        system("rm -Rf %s" % llhdir)

class GetLenLlh_V_Agg(Target):
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
        self.setFollowOnTarget(GetLenLlh_V_Cleanup(self.outfile))

class GetLenLlh_V(Target):
    def __init__(self, model, db_dir, v, l, outfile):
        Target.__init__(self)
        self.model = model
        self.db_dir = db_dir
        self.v = v
        self.l = l
        self.outfile = outfile
    
    def run(self):
        # split db clones by j
        js = []
        tempdir = "%s-db_jsplit" % self.outfile
        system("mkdir -p %s" % tempdir)
        for sam in os.listdir(self.db_dir):
            infile = os.path.join(self.db_dir, sam, self.v, self.l)
            if not os.path.exists(infile):
                continue
            clones = pickle.load(gzip.open(infile, "rb"))
            j2clones = lcommon.split_clones_by_j(clones)
            tempsamdir = os.path.join(tempdir, sam)
            system('mkdir -p %s' % tempsamdir)
            for j, jclones in j2clones.iteritems():
                tempjfile = os.path.join(tempsamdir, j)  #  l/v-db_jsplit/sam/jfile
                pickle.dump(jclones, gzip.open(tempjfile, 'wb'))
                if j not in js:
                    js.append(j)
        self.logToMaster("Done spliting clones by j for length %s V %s\n" %
                          (self.l, self.v))
        
        tempoutdir = "%s-llh_jsplit" % self.outfile  #  l/v-llh_jsplit/jllh
        system("mkdir -p %s" % tempoutdir)
        for j in js:
            joutfile = os.path.join(tempoutdir, j)
            self.addChildTarget(GetLenLlh_VJ(self.model, tempdir, j, joutfile))

        self.setFollowOnTarget(GetLenLlh_V_Agg(tempoutdir, self.outfile))

class GetLenLlhAgg(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        llh = 0.0
        for v in os.listdir(self.indir):
            vfile = os.path.join(self.indir, v)
            vllh = pickle.load(gzip.open(vfile, "rb"))
            llh += vllh
        if llh == 0:
            log_llh = float('-inf')
        else:
            log_llh = log10(llh)
        pickle.dump(log_llh, gzip.open(self.outfile, "wb"))

class GetLenLlh(Target):
    def __init__(self, db_dir, length, model, outdir):
        Target.__init__(self)
        self.db_dir = db_dir
        self.length = length
        self.model = model
        self.outdir = outdir
    
    def run(self):
        vs = get_union_children(self.db_dir)
        for v in vs:
            outfile = os.path.join(self.outdir, v)
            self.addChildTarget(GetLenLlh_V(self.model, self.db_dir, v,
                                            self.length, outfile))
        aggfile = os.path.join(self.outdir, "%s.pickle" % self.length)
        self.setFollowOnTarget(GetLenLlhAgg(self.outdir, aggfile))

class GetLenLlhsAgg(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
    
    def run(self):
        f = open(self.outfile, 'w')
        f.write("#Length\tLog_likelihood\n")
        lens = sorted([int(l) for l in os.listdir(self.indir)])
        for l in lens:
            lfile = os.path.join(self.indir, str(l), "%s.pickle" % str(l))
            llh = pickle.load(gzip.open(lfile, "rb"))
            f.write("%d\t%f\n" % (l, llh))
        f.close()

class GetLenLlhs(Target):
    '''compute the likelihood of all clones:
    global_dir/
            length1/
                v1
                v2
                ...
            length2/
            ...
        ...
    '''
    def __init__(self, db_dir, lens, model, outfile):
        Target.__init__(self)
        self.db_dir = db_dir
        self.lens = lens
        self.model = model
        self.outfile = outfile

    def run(self):
        self.logToMaster("Starting to compute llh for each length...\n")
        global_dir = self.getGlobalTempDir()
        for l in self.lens:
            outdir = os.path.join(global_dir, str(l))
            system("mkdir -p %s" % outdir)
            self.addChildTarget(GetLenLlh(self.db_dir, l, self.model, outdir))
        self.setFollowOnTarget(GetLenLlhsAgg(global_dir, self.outfile))

class Setup(Target):
    def __init__(self, db_dir, model_dir, outfile, options):
        Target.__init__(self)
        self.db_dir = db_dir
        self.model_dir = model_dir
        self.outfile = outfile
        self.options = options

    def run(self):
        model = rcommon.get_median_model(self.model_dir)
        self.logToMaster("Done computing median model.\n")
       
        lens = get_lens(self.db_dir)
        self.logToMaster("Done getting lengths: %s\n" %
                                             ",".join([str(l) for l in lens]))
        self.addChildTarget(GetLenLlhs(self.db_dir, lens, model, self.outfile))

def main():
    usage = "%prog <db_dir> <model_dir> <outfile> [options]"
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    db_dir = args[0]
    model_dir = args[1]
    outfile = args[2]

    i = Stack(Setup(db_dir, model_dir, outfile, options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.len_llh import *
    main()

