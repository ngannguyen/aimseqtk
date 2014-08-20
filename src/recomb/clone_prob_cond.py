#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Get the median models
Input a list of clones
Compute the probability of observing each clone
condition on length
Compute the number of expected samples in each group
that each clone is observed in
'''

import os
import sys
import numbers
import re
from math import log10
import cPickle as pickle
import gzip
import numpy as np
from scipy.stats import poisson

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon
import aimseqtk.src.recomb.recomb_common as rcommon


def read_clone_file(file, short=False):
    # if short=True, infile only has 2 fields: <clone>\t<numsam>
    clone2names = {}
    f = open(file, 'r')
    f.readline()
    for line in f:
        items = line.strip('\n').split('\t')
        #assert len(items) == 7
        clone = items[0]
        if short:
            clone2names[clone] = int(items[1])
        else:
            in_names = items[5].split(',')
            out_names = items[6].split(',')
            clone2names[clone] = (in_names, out_names)
    f.close()
    return clone2names

def get_ntclones(vaaj, sams, db_dir):
    # get all nt sequences that result in v_aa_j
    sam2ntclones = {}
    items = vaaj.split('_')
    assert len(items) == 3
    v = items[0]
    l = len(items[1])
    for sam in sams:
        if not sam:
            continue
        ntclones = []
        clonefile = os.path.join(db_dir, sam, v, str(l))
        clones = pickle.load(gzip.open(clonefile, "rb"))
        for clone in clones:
            if clone.get_vseqj() == vaaj:
                if clone.vdel is not None:
                    ntclones.append(clone)
        sam2ntclones[sam] = ntclones
    return sam2ntclones

def get_group_likelihood(log_p, allsams, num_present, sam2total):
    llhs = []
    for sam in allsams:
        logmu = log_p + log10(sam2total[sam])
        llh = poisson.logsf(1, 10 ** logmu)
        llhs.append(llh)
    sum_llh = sum([10 ** l for l in llhs])
    x = num_present
    return poisson.logsf(x, sum_llh), sum_llh

def get_lencount(db_dir):
    len2count = {}
    s = os.path.basename(db_dir)
    for v in os.listdir(db_dir):
        if v == s:
            continue
        vdir = os.path.join(db_dir, v)
        for l in os.listdir(vdir):
            lfile = os.path.join(vdir, l)
            clones = pickle.load(gzip.open(lfile, 'rb'))
            count = len(clones)
            l = int(l)
            if l not in len2count:
                len2count[l] = count
            else:
                len2count[l] += count
    return len2count

def get_lencount_fixedlen(lencount_dir, l):
    sam2count = {}
    for s in os.listdir(lencount_dir):
        sfile = os.path.join(lencount_dir, s)
        len2count = pickle.load(gzip.open(sfile, 'rb'))
        if l in len2count:
            sam2count[s] = len2count[l]
        else:
            sam2count[s] = 0
    return sam2count

def get_numclone_fixedlen(db_dir, l):
    # get number of clones with length l for each sample
    sam2numlen = {}
    for s in os.listdir(db_dir):
        sdir = os.path.join(db_dir, s)
        numlen = 0
        for v in os.listdir(sdir):
            if v == s:
                continue
            infile = os.path.join(db_dir, s, v, str(l))
            if os.path.exists(infile):
                clones = pickle.load(gzip.open(infile, 'rb'))
                numlen += len(clones)
        sam2numlen[s] = numlen
    return sam2numlen

def aaclone_llh(clone, cloneinfo, model, lencount_dir, group2sams, 
                outfile, ingroup, outgroup, len2llh, aa_llh):
    f = open(outfile, 'w')
    f.write("sample\tnum_ntclones\tprob_observed\n")
        
    items = clone.split('_')
    v = items[0]
    l = len(items[1])

    len_llh = len2llh[l]
    aa_llh_cond = aa_llh - len_llh
    f.write("#Log_len_llh: %f\n" % len_llh)
    f.write("#Log_aaclone_llh_cond: %f\n" % aa_llh_cond)
    #sam2numlen = get_numclone_fixedlen(db_dir, l)
    sam2numlen = get_lencount_fixedlen(lencount_dir, l)

    #for i, sams in enumerate([insams, outsams]):
    #    sam2ntclones = get_ntclones(clone, sams, db_dir)
    #    f.write("#Group_%d\n" % (i + 1))
    #    for sam, ntclones in sam2ntclones.iteritems():
    #        totallen = sam2numlen[sam]
    #        f.write("%s\t%d\n" % (sam, len(ntclones)))
    
    # calc prob to observe the aa clones (sum of all nt events)

    avr_totallen = sum(sam2numlen.values()) / len(sam2numlen)

    avr_logmu = aa_llh_cond + log10(avr_totallen)
    avr_aa_llh = poisson.logsf(1, 10 ** avr_logmu)
    f.write("#Clone_log_likelihood: %f, %f\n" % (aa_llh_cond, avr_aa_llh))
    
    if isinstance(cloneinfo, int):
        obs_numsam = cloneinfo  # observed # samples
        allsams = []
        for sams in group2sams.values():
            allsams.extend(sams)
        group_llh, expected_sams = get_group_likelihood(aa_llh_cond, allsams,
                                                        obs_numsam, sam2numlen)
        f.write("#Llh,Obs vs Exp:\t%f\t%d\t%f\n#\n" % (group_llh, obs_numsam,
                                                       expected_sams))
    else:
        insams = cloneinfo[0]
        outsams = cloneinfo[1]
        ingroup_llh, in_expected_sams = get_group_likelihood(aa_llh_cond,
                                                        group2sams[ingroup],
                                                        len(insams), sam2numlen)
        outgroup_llh, out_expected_sams = get_group_likelihood(aa_llh_cond,
                                                      group2sams[outgroup],
                                                      len(outsams), sam2numlen)
        f.write("#Ingrp vs Outgrp: %f vs %f\n" % (ingroup_llh, outgroup_llh))
        f.write("#Expected Ingrp, Outgrp:\t%f\t%f\n" % (in_expected_sams,
                                                           out_expected_sams))
        f.write("#Observed Ingrp, Outgrp:\t%d\t%d\n#\n" %(len(insams),
                                                          len(outsams)))
    f.close()

def read_clonesize(file):
    s2clones = {}
    g2sams = {}
    f = open(file, 'r')
    f.readline()
    
    sams = []
    for line in f:
        items = line.strip().split('\t')
        sample = items[0]
        clones = int(float(items[1]))
        if re.search("Avr", sample):
            group = sample.replace("_Avr", "")
            g2sams[group] = sams
            sams = []
            s2clones[group] = clones
        else:
            s2clones[sample] = clones
            sams.append(sample)
    f.close()
    return s2clones, g2sams

def read_llh(file, intkey=False):
    len2llh = {}
    f = open(file, 'r')
    for line in f:
        if not line or line[0] == '#':
            continue
        items = line.strip().split('\t')
        l = items[0]
        if intkey:
            l = int(items[0])
        llh = float(items[1])
        len2llh[l] = llh
    f.close()
    return len2llh

class AaCloneLlh(Target):
    def __init__(self, clone, cloneinfo, model, lencount_dir, group2sams,
                 outfile, ingroup, outgroup, len2llh, clonellh):
        Target.__init__(self)
        self.clone = clone
        self.cloneinfo = cloneinfo
        #self.insams = sams[0]
        #self.outsams = sams[1]
        self.model = model
        self.lencount_dir = lencount_dir
        self.group2sams = group2sams
        self.outfile = outfile
        self.ingroup = ingroup
        self.outgroup = outgroup
        self.len2llh = len2llh
        self.clonellh = clonellh

    def run(self):
        #self.logToMaster("AaCloneLllh: %s\n" % self.clone)
        aaclone_llh(self.clone, self.cloneinfo, self.model,
                    self.lencount_dir, self.group2sams, self.outfile,
                    self.ingroup, self.outgroup, self.len2llh, self.clonellh)
        #self.logToMaster("DONE AaCloneLllh: %s\n" % self.clone)

class Setup(Target):
    def __init__(self, args):
        Target.__init__(self)
        self.clone_file = args[0]
        self.model = args[1]
        self.db_dir = args[2]
        self.numclone_file = args[3]
        self.outdir = args[4]
        self.ingroup = args[5]
        self.outgroup = args[6]
        self.lenllh = args[7]
        self.clonellh = args[8]
    
    def run(self):
        system("mkdir -p %s" % self.outdir)
        clone2sams = read_clone_file(self.clone_file, True)
        if os.path.isdir(self.model):
            model = rcommon.get_median_model(self.model)
        else:
            model = pickle.load(gzip.open(self.model, 'rb'))
        sam2total, group2sams = read_clonesize(self.numclone_file)
        len2llh = read_llh(self.lenllh, intkey=True)
        clone2llh = read_llh(self.clonellh)
        
        global_dir = self.getGlobalTempDir()
        lencount_dir = os.path.join(global_dir, "sam2len2count")
        system("mkdir -p %s" % lencount_dir)
        for s in os.listdir(self.db_dir):
            samdir = os.path.join(self.db_dir, s)
            lencount_file = os.path.join(lencount_dir, s)
            self.addChildTarget(GetLencount(samdir, lencount_file))
        self.setFollowOnTarget(GetLlhs(clone2sams, self.outdir, model,
                                       lencount_dir, group2sams, self.ingroup,
                                       self.outgroup, len2llh, clone2llh))

class GetLlhs(Target):
    def __init__(self, clone2sams, outdir, model, lencount_dir, group2sams,
                 ingroup, outgroup, len2llh, clone2llh):
        Target.__init__(self)
        self.clone2sams = clone2sams
        self.outdir = outdir
        self.model = model
        self.lencount_dir = lencount_dir
        self.group2sams = group2sams
        self.ingroup = ingroup
        self.outgroup = outgroup
        self.len2llh = len2llh
        self.clone2llh = clone2llh
    
    def run(self):
        for clone, sams in self.clone2sams.iteritems():
            outfile = os.path.join(self.outdir, clone)
            clonellh = self.clone2llh[clone]
            self.addChildTarget(AaCloneLlh(clone, sams, self.model,
                                   self.lencount_dir,
                                   self.group2sams, outfile, self.ingroup,
                                   self.outgroup, self.len2llh, clonellh))

class GetLencount(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        l2c = get_lencount(self.indir)
        pickle.dump(l2c, gzip.open(self.outfile, 'wb'))

def main():
    usage = ("%prog <clone_file> <model_dir|model_pickle_file> <db_dir> <numclone_file>" +
             "<out_dir> <ingroup> <outgroup> <lenllh> <clonellh>")
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    
    i = Stack(Setup(args)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)
    
if __name__ == '__main__':
    from aimseqtk.src.recomb.clone_prob_cond import *
    main()

