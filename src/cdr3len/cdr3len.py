#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''CDR3 length distribution: Reads and Clones
% reads/% clones (yaxis) vs CDR3 aa_length (xaxis) 
ttests: 1/ for the median length of group1 vs group2
        2/ for each length, compare the proportion (reads/clones)
        in group1 vs group2 having that length
'''

import os
import sys
import gzip
import cPickle as pickle
from optparse import OptionGroup

from numpy import median
from jobTree.scriptTree.target import Target
from sonLib.bioio import system

from aimseqtk.lib.statcommon import SampleStat
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
from aimseqtk.lib.common import StatAnalyses
from aimseqtk.lib.common import Analysis
import aimseqtk.lib.statcommon as statcommon
import aimseqtk.src.cdr3len.cdr3len_plot as ldplot


class LenDistStat(SampleStat):
    def __init__(self):
        SampleStat.__init__(self)
        self.len2clones = {}
        self.len2reads = {}
        self.median_clones = None  # median cdr3 length based on # clones
        self.median_reads = None  # median cdr3 length based on # reads

    def set_stats(self, len2clones, len2reads):
        self.__setitem__('len2clones', len2clones)
        clonevec = statcommon.sizedict_to_vec(self.len2clones, True)
        self.__setitem__('median_clones', median(clonevec))
        self.len2reads = len2reads
        readvec = statcommon.sizedict_to_vec(self.len2reads, True)
        self.median_reads = median(readvec)

def sample_lendist_stat(sample, samdir, args=None):
    # lendist with counts and with number of clones
    len2clones = {}
    len2reads = {}
    clones = libsample.sample_all_clones(samdir)
    totalclone = 0
    for clone in clones:
        if clone.aa:
            l = len(clone.aa)
            if clone.vdel is not None:
                totalclone += 1
                if l not in len2clones:
                    len2clones[l] = 1
                else:
                    len2clones[l] += 1
            if l not in len2reads:
                len2reads[l] = clone.freq
            else:
                len2reads[l] += clone.freq
    # convert the number of clones into % total clones
    for l, numclone in len2clones.iteritems():
        len2clones[l] = float(numclone) / totalclone

    stat = LenDistStat()
    stat.set_sample_info(sample)
    stat.set_stats(len2clones, len2reads)
    return stat

def get_lens(stats):
    # return the union list of cdr3 lengths of all samples
    lens = []
    for stat in stats:
        ls = stat.len2clones.keys()
        for l in ls:
            if l not in lens:
                lens.append(l)
    return sorted(lens)

def lendist_median(name2obj, attr, outfile):
    l2freqs = {}
    lens = get_lens(name2obj.values())
    for l in lens:
        l2freqs[l] = []
    for obj in name2obj.values():
        attrstr = "len2" + attr
        l2v = obj[attrstr]
        for l in lens:
            if l not in l2v:
                l2freqs[l].append(0.0)
            else:
                l2freqs[l].append(l2v[l])
    l2medfreq = {}
    for l, freqs in l2freqs.iteritems():
        l2medfreq[l] = median(freqs)
    
    # print to outfile
    f = open(outfile, 'w')
    f.write("#Length\tFrequency\n")
    for l in sorted(lens):
        f.write("%d\t%f\n" % (l, l2medfreq[l]))
    f.close()

def lendist_ttests(attr, outfile, g2n, name2obj, matched, pcutoff):
    testtype = statcommon.get_test_type(g2n, matched)
    f = open(outfile, 'w')
    # compare median length
    f.write(("#Category\tGroup1_Group2\tt_val\tp_val\tMean1 +/- Std1\t" +
             "Mean2 +/- Std2\n"))
    medattr = "median_" + attr 
    pair2tp, group2mean = statcommon.ttest_allpairs(g2n, name2obj, matched,
                                                    medattr, testtype=testtype)
    statcommon.ttest_write(f, medattr, pair2tp, group2mean)

    # compare freq of each length
    f.write("#Length\n")
    lens = get_lens(name2obj.values())
    lenattr = "len2" + attr
    for l in lens:
        pair2tp, group2mean = statcommon.ttest_allpairs(g2n, name2obj,
                                  matched, attr=None,
                                  func=statcommon.obj_dictattr_lookup,
                                  func_args=(lenattr, l), testtype=testtype)
        statcommon.ttest_write(f, str(l), pair2tp, group2mean, pcutoff)
    f.close()

class LenDistAnalyses(StatAnalyses):
    '''Draw lendist plots and perform ttests
    '''
    def __init__(self, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)

    def run(self):
        #name2obj = libcommon.load_pickledir_to_dict(self.indir)
        self.load_indir()
        name2obj = self.name2obj
        attrs = ['clones', 'reads'] 
        plotfmt = self.opts.plotformat
        for attr in attrs:
            plotfile = os.path.join(self.outdir, "%s" % attr)
            if len(name2obj) >= 20:
                ldplot.draw_lendist_avr(name2obj, attr, plotfile, plotfmt,
                                    self.opts.dpi)
            else:
                ldplot.draw_lendist(name2obj, attr, plotfile, plotfmt,
                                    self.opts.dpi)
            # ttests
            g2n = self.opts.group2samples
            if g2n:
                ttestfile = os.path.join(self.outdir, "ttests_%s.txt" % attr)
                lendist_ttests(attr, ttestfile, g2n, name2obj,
                               self.opts.matched, self.opts.pval)

            # median len distribution
            meddist_file  = os.path.join(self.outdir, "median_lendist_%s.txt" % attr)
            lendist_median(name2obj, attr, meddist_file)

class LenDist(Analysis):
    '''Set up children jobs to compute CDR3 length distribution for
    each sample, and set follow on target to make plots and do ttests
    '''
    def __init__(self, indir, outdir, opts):
        Analysis.__init__(self, indir, outdir, opts)

    def run(self):
        opts = self.opts
        global_dir = self.getGlobalTempDir()
        ld_dir = os.path.join(global_dir, "lendist_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % ld_dir)
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            sample = pickle.load(gzip.open(os.path.join(samdir, sam), 'rb'))
            outfile = os.path.join(ld_dir, "%s.pickle" % sam)
            self.addChildTarget(libsample.SampleAnalysis(sample, samdir,
                                                 outfile, sample_lendist_stat))
        self.setFollowOnTarget(LenDistAnalyses(ld_dir, self.outdir, opts))

