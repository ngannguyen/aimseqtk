#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
For matched samples:
1/ Track a specific clones over different time points
Or
2/ Track clones >= minFreq:
    Get all clones with freq >= minFreq for each sample at each time point
    Track frequency of those clones across all time points
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
from aimseqtk.lib.common import StatAnalyses
from aimseqtk.lib.common import Analysis
import aimseqtk.lib.statcommon as statcommon
import aimseqtk.src.overlap.trackclone_plot as tcplot


def check_track_clone_options(parser, options):
    options.clones = None
    if options.clonefile:
        libcommon.check_options_file(options.clonefile)
        options.clones = libcommon.read_list(options.clonefile)
    if options.track_minfreq > 1.0:
        parser.error("--track_minfreq is the proportion, not %, should be < 1")

def add_track_clone_options(parser):
    group = OptionGroup(parser, "Clone-tracking options")
    group.add_option('--track_minfreq', dest='track_minfreq', default=0.1,
                     help=('Clones with freqs >= track_minfreq are tracked.'))
    group.add_option('--track_clones', dest='clonefile',
                     help=('File that lists clones for tracking'))
    parser.add_option_group(group)
    
def get_matched_samples(sample, group, groups, group2samples):
    # return the matched samples of "sample" sorted by "groups"
    index = -1
    matched_samples = []
    for i, s in enumerate(group2samples[group]):
        if s == sample:
            index = i
            break
    if index < 0:
        raise ValueError("Group %s does not have Sample %s" % (group, sample))
    for g in groups:
        if g not in group2samples:
            raise ValueError("Group %s is not in group2samples" % group)
        if len(group2samples[g]) <= index:
            raise ValueError("Group %s does not have enough members" % group)
        matched_samples.append(group2samples[g][index])
    return matched_samples, index

def clone_get_samples(cloneid, samples, sizetype='count'):
    # get all samples that contain "clone"
    name2size = {}
    for sample in samples:
        for clone in sample.clones:
            ids = clone.get_vjseq_ids()
            for id in ids:
                if id == cloneid:
                    size = clone[sizetype] / len(ids)
                    if size == 0 and sizetype == 'count':
                        size = 1
                    if sample.name not in name2size:
                        name2size[sample.name] = size
                    else:
                        name2size[sample.name] += size
                    break
    return name2size

def track_clone(clone, sample2size, groups, group2samples, sample2group):
    # tracking clone's frequences over different "groups" in the group
    # order given. sample2size are all samples containing "clone"
    indices = []
    rows = []  # each row = each matched set of samples
    visited = []  #samples already visited
    samples = sorted(sample2size.keys())
    for sample in samples:
        if sample in visited:
            continue
        group = sample2group[sample]
        matched_samples, index = get_matched_samples(sample, group, groups,
                                                                 group2samples)
        row = []
        for s in matched_samples:
            visited.append(s)
            if s in sample2size:
                row.append(sample2size[s])
            else:
                row.append(0.0)
        rows.append(row)
        indices.append(index)
    return rows, indices

def sample_top_clones(sample, minsize, attr='freq'):
    topclone2size = {}
    for clone in sample.clones:
        ids = clone.get_vjseq_ids()
        size = clone[attr] / len(ids)
        if size == 0 and attr == 'count':
            size = 1
        if size >= minsize:
            for id in ids:
                if id in topclone2size:
                    topclone2size[id] += size
                else:
                    topclone2size[id] = size
    return topclone2size

def top_clones(samples, minsize, attr='freq'):
    topclone2name2size = {}
    for sample in samples:
        topclone2size = sample_top_clones(sample, minsize, attr)
        for tc, size in topclone2size.iteritems():
            if tc not in topclone2name2size:
                topclone2name2size[tc] = {sample.name: size}
            else:
                topclone2name2size[tc][sample.name] = size
    return topclone2name2size

def track_top_clones(samples, minsize, groups, group2samples, sample2group,
                                                                  attr='freq'):
    tc2rows = {}
    tc2indices = {}
    tc2name2size = top_clones(samples, minsize, attr)
    for tc in tc2name2size:
        name2size = clone_get_samples(tc, samples, attr) 
        rows, indices = track_clone(tc, name2size, groups, group2samples, sample2group)
        tc2rows[tc] = rows
        tc2indices[tc] = indices
    return tc2rows, tc2indices

def tab_track_clones(clone2rows, groups, outfile):
    f = open(outfile, 'w')
    f.write("#Clone\tIndex\t%s\n" % "\t".join(groups))
    for clone, (rows, indices) in clone2rows.iteritems():
        for i, row in enumerate(rows):
            index = indices[i]
            rowstr = "\t".join(["%.2e" % c for c in row])
            f.write("%s\t%d\t%s\n" % (clone, index, rowstr))
    f.close()

#====== JobTree Targets ======
class TrackClonesAgg(StatAnalyses):
    '''
    '''
    def __init__(self, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)
    
    def run(self):
        n2o = self.name2obj
        outfile = os.path.join(self.outdir, "trackclones.txt")
        tab_track_clones(n2o, opts.groups, outfile)
        if self.opts.makeplots:
            plotdir = os.path.join(self.outdir, "plots")
            system('mkdir -p %s' % plotdir)
            for clone, (rows, indices) in n2o.iteritems():
                outbase = os.path.join(plotdir, clone)
                tcplot.draw_track_clone(clone, rows, self.opts.groups,
                                        outbase, self.opts)

class TrackClone(Target):
    '''track a specific clone
    '''
    def __init__(self, clone, samples, outfile, opts, s2g):
        Target.__init__(self)
        self.clone = clone
        self.samples = samples
        self.outfile = outfile
        self.opts = opts
        self.s2g = s2g

    def run(self):
        name2size = clone_get_samples(self.clone, self.samples, attr='freq')
        rows, indices = track_clone(self.clone, name2size, self.opts.groups,
                                    self.opts.group2samples, g2s, self.s2g)
        pickle.dump((rows, indices), gzip.open(self.outfile, "wb"))

class TrackClones(Analysis):
    '''Track a number of input clones
    '''
    def __init__(self, samples, outdir, opts, clones):
        Analysis.__init__(self, samples, outdir, opts)
        self.clones = clones
    
    def run(self):
        opts = self.opts
        sample2group = libcommon.get_val2key_1to1(opts.group2samples)
        global_dir = self.getGlobalTempDir()
        for clone in self.clones:
            outfile = os.path.join(global_dir, "%s.pickle" % clone)
            self.addChildTarget(TrackClone(clone, self.samples, outfile, opts,
                                           sample2group))
        self.setFollowOnTarget(TrackClonesAgg(global_dir, self.outdir, opts))

class TrackTopClones(TrackClones):
    '''Set up children jobs to track top clones
    '''
    def __init__(self, samples, outdir, opts):
        topclones = top_clones(samples, opts.track_minfreq)
        Analysis.__init__(self, samples, outdir, opts, topclones)

