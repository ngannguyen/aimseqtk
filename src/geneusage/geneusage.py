#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''Gene Usage Distribution
V, J, D, VJ, VDJ
ttests
'''

import os
import sys
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

from aimseqtk.lib.statcommon import SampleStat
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
from aimseqtk.lib.common import Analysis
from aimseqtk.lib.common import StatAnalyses
import aimseqtk.lib.statcommon as statcommon
import aimseqtk.src.geneusage.geneusage_plot as guplot


class GeneUsageStat(SampleStat):
    '''Contain geneusage stats of a specific sample
    '''
    def __init__(self):
        SampleStat.__init__(self)
        self.type2gene2clones = {}
        self.type2gene2reads = {}
        types = ['v', 'd', 'j', 'vj', 'dj']
        for type in types:
            self.type2gene2clones[type] = {}
            self.type2gene2reads[type] = {}

    def set_stats(self, type2gene2clones, type2gene2reads):
        self.type2gene2clones = type2gene2clones
        self.type2gene2reads = type2gene2reads

def sample_geneusage_stat(sample, samdir, args=None):
    # gene usage of a specific number
    # initialize usage
    types = ['v', 'd', 'j', 'vj', 'dj']
    type2gene2clones = {}
    type2gene2reads = {}
    for t in types:
        type2gene2clones[t] = {}
        type2gene2reads[t] = {}

    # get usage
    clones = libsample.sample_all_clones(samdir)
    for clone in clones:
        # update each genetype usage
        genetypes = ['v', 'd', 'j']
        for type in genetypes:
            gene2clones = type2gene2clones[type]
            gene2reads = type2gene2reads[type]
            gene = clone[type]
            numclone = 1.0
            freq = clone.freq
            if gene not in gene2clones:
                gene2clones[gene] = numclone
                gene2reads[gene] = freq
            else:
                gene2clones[gene] += numclone
                gene2reads[gene] += freq

        # update gene combination usage
        for type in ['vj', 'dj']:
            gene2clones = type2gene2clones[type]
            gene2reads = type2gene2reads[type]
            g0 = clone[type[0]]
            g1 = clone[type[1]]
            numclone = 1.0
            freq = clone.freq
            combi = "|".join([g0, g1])
            if combi not in gene2clones:
                gene2clones[combi] = numclone
                gene2reads[combi] = freq
            else:
                gene2clones[combi] += numclone
                gene2reads[combi] += freq
    # convert the number of clones into % total clones:
    for type, gene2clones in type2gene2clones.iteritems():
        for gene, numclone in gene2clones.iteritems():
            gene2clones[gene] = float(numclone) / sample.numclone

    # get the stat obj
    stat = GeneUsageStat()
    stat.set_sample_info(sample)
    stat.set_stats(type2gene2clones, type2gene2reads)
    return stat

def get_genes(stats, genetype):
    # return the union list of genes of all samples
    genes = []
    for stat in stats:
        currgenes = stat.type2gene2clones[genetype].keys()
        for g in currgenes:
            if g not in genes:
                genes.append(g)
    return genes

def get_geneusage(stat, args):
    # return sample gene usage of a specific category
    assert len(args) == 3
    count_type = args[0]
    gene_type = args[1]
    gene = args[2]
    attr = "type2gene2" + count_type 
    assert attr in stat.getitems()
    assert gene_type in stat[attr]
    if gene not in stat[attr][gene_type]:
        return 0.0
    return stat[attr][gene_type][gene]

def geneusage_ttests(attr, type, outfile, g2n, name2obj, matched, pcutoff):
    # attr: "clones"/"reads"; type=genetype: v, d, j, vj, dj
    f = open(outfile, 'w')
    f.write(("Name\tGroup1_Group2\tt_val\tp_val\tMean1 +/- Std1\t" +
             "Mean2 +/- Std2\n"))
    genes = get_genes(name2obj.values(), type)
    for gene in genes:
        pair2tp, group2mean = statcommon.ttest_allpairs(g2n, name2obj, matched,
                                                 attr=None, func=get_geneusage,
                                                 func_args=(attr, type, gene))
        statcommon.ttest_write(f, gene, pair2tp, group2mean, pcutoff)
    f.close()

class GeneUsageTtests(Target):
    '''Perform ttests for a specific gene category (v, d, j, vj or dj)
    for each pair of groups
    '''
    def __init__(self, attr, type, outfile, g2n, name2obj, matched, pcutoff):
        Target.__init__(self)
        self.attr = attr
        self.type = type
        self.outfile = outfile
        self.g2n = g2n
        self.name2obj = name2obj
        self.matched = matched
        self.pcutoff = pcutoff

    def run(self):
        geneusage_ttests(self.attr, self.type, self.outfile, self.g2n,
                         self.name2obj, self.matched, self.pcutoff)

class GeneUsageAnalyses(StatAnalyses):
    '''Draw usage plots
    Perform ttests
    '''
    def __init__(self, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)
    
    def run(self):
        self.load_indir()
        attrs = ['clones', 'reads']
        types = ['v', 'd', 'j', 'vj', 'dj']
        
        g2n = self.opts.group2samples
        for attr in attrs:
            # draw usage plots
            if self.opts.makeplots:
                plotdir = os.path.join(self.outdir, "%s_plots" % attr)
                system("mkdir -p %s" % plotdir)
                for type in types:
                    plotfile = os.path.join(plotdir, type)
                    genes = get_genes(self.name2obj.values(), type)
                    self.addChildTarget(guplot.GeneUsagePlot(self.name2obj,
                                       attr, type, genes, plotfile, self.opts))
            # ttests
            if g2n:
                ttestdir = os.path.join(self.outdir, "%s_ttests" % attr)
                system("mkdir -p %s" % ttestdir)
                for type in types:
                    ttestfile = os.path.join(ttestdir, "%s.txt" % type)
                    self.addChildTarget(GeneUsageTtests(attr, type, ttestfile,
                       g2n, self.name2obj, self.opts.matched, self.opts.pval))

class GeneUsage(Analysis):
    '''Set up children jobs to compute gene usage for each sample
    and follow on job to do geneusage analyses, including makeing
    plots and performing ttests
    '''
    def __init__(self, indir, outdir, opts):
        Analysis.__init__(self, indir, outdir, opts)

    def run(self):
        global_dir = self.getGlobalTempDir()
        gu_dir = os.path.join(global_dir, "geneusage_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % gu_dir)
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            sample = pickle.load(gzip.open(os.path.join(samdir, sam), 'rb'))
            outfile = os.path.join(gu_dir, "%s.pickle" % sam)
            self.addChildTarget(libsample.SampleAnalysis(sample, samdir,
                                            outfile, sample_geneusage_stat))
        self.setFollowOnTarget(GeneUsageAnalyses(gu_dir, self.outdir,
                                                                self.opts))




