#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Get the median models
Input a list of clones
Compute the probability of observing each clone
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
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon
import aimseqtk.src.recomb.recomb_common as rcommon


def read_clone_file(file):
    clone2names = {}
    f = open(file, 'r')
    f.readline()
    for line in f:
        items = line.strip('\n').split('\t')
        assert len(items) == 7
        clone = items[0]
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
    for sam in sams:
        if not sam:
            continue
        ntclones = []
        clonefile = os.path.join(db_dir, sam, v)
        clones = pickle.load(gzip.open(clonefile, "rb"))
        for clone in clones:
            if clone.get_vseqj() == vaaj:
                if clone.vdel is not None:
                    ntclones.append(clone)
        sam2ntclones[sam] = ntclones
    return sam2ntclones

def n_choose_k(n, k):
    numerator = factorial(n)
    denominator = (factorial(k) * factorial(n - k))
    return numerator / denominator

def get_group_likelihood(log_p, allsams, presentsams, sam2total):
    llhs = []
    for sam in allsams:
        logmu = log_p + log10(sam2total[sam])
        llh = poisson.logsf(1, 10 ** logmu)
        llhs.append(llh)
    sum_llh = sum([10 ** l for l in llhs])
    x = len(presentsams)
    return poisson.logsf(x, sum_llh)

def aaclones_likelihood(clone2sams, model, db_dir, sam2total, group2sams,
                        outfile, ingroup, outgroup):
    f = open(outfile, 'w')
    f.write("sample\tnum_ntclones\tprob_observed\n")
    for clone, (insams, outsams) in clone2sams.iteritems():
        f.write("#%s\n" % clone)
        events = []
        event_llhs = []
        for i, sams in enumerate([insams, outsams]):
            if not sams:
                continue
            sam2ntclones = get_ntclones(clone, sams, db_dir)
            f.write("#Group_%d\n" % (i + 1))
            for sam, ntclones in sam2ntclones.iteritems():
                total = sam2total[sam]
                llhoods = []
                for ntclone in ntclones:
                    clonellhood = rcommon.ntclone_likelihood(ntclone, model)
                    #prob_observed = clonellhood + log10(total)
                    logmu = log10(total) + clonellhood
                    prob_observed = poisson.logsf(1, 10 ** logmu)  # prob. observing >=1 ntclone
                    llhoods.append(prob_observed)

                    if not rcommon.visited_event(events, ntclone):
                        events.append(ntclone)
                        event_llhs.append(clonellhood)
                        #if clonellhood != float(-inf):
                        #    event_llhs.append(clonellhood)

                llhoods_str = ",".join(["%f" % llh for llh in llhoods])
                f.write("%s\t%d\t%s\n" % (sam, len(ntclones), llhoods_str))
        
        # calc prob to observe the aa clones (sum of all nt events)
        if sum([10**llh for llh in event_llhs]) > 0:
            aa_llh = log10(sum([10**llh for llh in event_llhs]))
            avr_total = (sam2total[ingroup] + sam2total[outgroup]) / 2
            avr_logmu = aa_llh + log10(avr_total)
            avr_aa_llh = poisson.logsf(1, 10 ** avr_logmu)
            f.write("#Clone_log_likelihood: %f, %f\n" % (aa_llh, avr_aa_llh))
            
            ingroup_llh = get_group_likelihood(aa_llh, group2sams[ingroup],
                                               insams, sam2total)
            outgroup_llh = get_group_likelihood(aa_llh, group2sams[outgroup],
                                                outsams, sam2total)
            f.write("#Ingrp vs Outgrp: %f vs %f\n#\n" % (ingroup_llh, outgroup_llh))

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

def main():
    clone_file = sys.argv[1]
    model_dir = sys.argv[2]
    db_dir = sys.argv[3]
    numclone_file = sys.argv[4]
    outfile = sys.argv[5]
    ingroup = sys.argv[6]
    outgroup = sys.argv[7]

    clone2sams = read_clone_file(clone_file)
    model = rcommon.get_median_model(model_dir)
    sam2total, group2sams = read_clonesize(numclone_file)
    aaclones_likelihood(clone2sams, model, db_dir, sam2total, group2sams,
                        outfile, ingroup, outgroup)

if __name__ == '__main__':
    main()

