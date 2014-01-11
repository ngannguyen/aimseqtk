#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Clonesize distribution: discrete and cumulative
% total clones vs % total reads (relative clone size)
% total reads (reads of all clones) vs % total reads (clone size)
% total reads vs clone rank
Clones whose freq >= minFreq
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
import aimseqtk.lib.tabcommon as tabcommon
import aimseqtk.src.properties.clonesize_plot as csplot


class CloneSizeStat(SampleStat):
    def __init__(self, freqs=[]):
        SampleStat.__init__(self)
        if freqs:
            self.freqs = freqs
        else:
            #self.freqs = [0, 0.001, 0.01, 0.1, 1]
            self.freqs = [0, 0.00001, 0.0001, 0.001, 0.01]
        self.numclones = [0.0] * len(self.freqs)
        self.counts = [0.0] * len(self.freqs)
        self.topfreqs = []  # frequencies of the top clones
        self.numclones_cumul = self.numclones  # cumulative
        self.counts_cumul = self.counts  # cumulative
        self.topfreqs_cumul = []  # cumulative freqs of top clones

def get_attr_colfields(attr, obj):
    if attr in ['numclones', 'counts', 'numclones_cumul', 'counts_cumul']:
        return [str(f) for f in obj.freqs]
    else:
        return [str(i + 1) for i in range(len(obj.topfreqs))]

def sample_clonesize_stat(sample, samdir, freqs=[], numtop=50, args=None):
    # calculate: number of clones/ counts that lie within each freq
    # range. Note: freqs must be sorted, or the func will sort it
    # <numtop>: number of top clones whose freqs will be report
    stat = CloneSizeStat(sorted(freqs))
    stat.set_sample_info(sample)
    if args:
        numtop = args[0]
    clones = libsample.sample_all_clones(samdir)
    sorted_clones = sorted(clones, reverse=True, key=lambda c: c.freq)
   
    for index, clone in enumerate(sorted_clones):
        if index < numtop:
            stat.topfreqs.append(clone.freq)

        for i, minfreq in enumerate(stat.freqs):
            maxfreq = float('inf')
            if i + 1 < len(stat.freqs):
                maxfreq = stat.freqs[i + 1]
            if minfreq <= clone.freq and clone.freq < maxfreq:
                stat.numclones[i] += 1
                stat.counts[i] += clone.count
    # convert to frequencies:
    stat.numclones = [libcommon.get_pc(c, stat.numclone) for c in
                                                                stat.numclones]
    stat.counts = [libcommon.get_pc(c, stat.size) for c in stat.counts]
    # get cumulative stats:
    stat.numclones_cumul = libcommon.get_cumulative(stat.numclones)
    stat.counts_cumul = libcommon.get_cumulative(stat.counts)
    stat.topfreqs_cumul = libcommon.get_cumulative(stat.topfreqs, True)
    return stat

def add_clonesize_options(parser):
    group = OptionGroup(parser, "Clone size distribution analyses options")
    group.add_option('--cs_topclone', type='int', default=50,
                     help=('Number of top clones to report frequencies. ' +
                           'Default=%default'))
    group.add_option('--cs_force_all', action='store_true', default=False,
                     help=('By default, if the number of samples > 100, ' +
                           'draw group averages only. If specified, will ' +
                           'force to draw all samples.'))
    parser.add_option_group(group)


class CloneSizePlots(StatAnalyses):
    '''Make tables and plots
    '''
    def __init__(self, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)

    def run(self):
        self.load_indir()
        name2obj = self.name2obj
        opts = self.opts
        assert len(name2obj) > 0
        obj0 = name2obj.values()[0]
        numsam = len(name2obj)
        g2s = opts.group2samples
        g2name_avr = None
        g2avr = None
        if g2s:
            g2name_avr = libcommon.get_group_avr(name2obj, g2s)
            if opts.makeplots:
                g2avr = {}
                for g, na in g2name_avr.iteritems():
                    groupname = na[0]
                    groupavr = na[1]
                    g2avr[groupname] = groupavr
            
        attrs = ['numclones', 'counts', 'topfreqs', 'numclones_cumul',
                 'counts_cumul', 'topfreqs_cumul']
        txtdir = os.path.join(self.outdir, "txt_tables")
        system("mkdir -p %s" % txtdir)
        plotdir = None
        if opts.makeplots:
            plotdir = os.path.join(self.outdir, "figures")
            system("mkdir -p %s" % plotdir)
        
        for attr in attrs:
            txtfile = os.path.join(txtdir, "%s.txt" % attr)
            colfields = get_attr_colfields(attr, obj0)
            tabcommon.table(name2obj, txtfile, colfields, g2name_avr,
                                               g2s, keyattr=attr, islist=True)
            if opts.makeplots:
                plotfile = os.path.join(plotdir, attr)
                if numsam < 100 or opts.cs_force_all:
                    csplot.draw_clonesize_dist(name2obj, attr, plotfile,
                                                    opts.plotformat, opts.dpi)
                elif g2avr:
                    csplot.draw_clonesize_dist(g2avr, attr, plotfile,
                                                        opts.plotformat, opts.dpi)

class CloneSize(Analysis):
    '''Make plots of different clonesize distributions
    '''
    def __init__(self, indir, outdir, opts):
        Analysis.__init__(self, indir, outdir, opts)

    def run(self):
        opts = self.opts
        topclone = opts.cs_topclone
        global_dir = self.getGlobalTempDir()
        cs_dir = os.path.join(global_dir, "clonesize_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % cs_dir)
        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            sample = pickle.load(gzip.open(os.path.join(samdir, sam), 'rb'))
            outfile = os.path.join(cs_dir, "%s.pickle" % sam)
            self.addChildTarget(libsample.SampleAnalysis(sample, samdir,
                                     outfile, sample_clonesize_stat, topclone))
        self.setFollowOnTarget(CloneSizePlots(cs_dir, self.outdir, opts))

