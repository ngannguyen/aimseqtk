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
                     type='float',
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

def clone_get_samples(cloneid, indir, vj, sizetype='count'):
    # get all samples that contain "clone"
    name2size = {}
    for sam in os.listdir(indir):
        vjfile = os.path.join(indir, sam, vj)
        if os.path.exists(vjfile):
            clones = pickle.load(gzip.open(vjfile, 'rb'))
            for clone in clones:
                id = clone.get_vseqj()
                if id == cloneid:
                    if sam not in name2size:
                        name2size[sam] = clone[sizetype]
                    else:
                        name2size[sam] += clone[sizetype]
                    break
    return name2size

def track_clone_no_matched(clone, sample2size, groups, group2samples):
    rows = []
    samples = []
    for group in groups:
        row = []
        for sample in group2samples[group]:
            if sample in sample2size:
                samples.append(sample)
                row.append(sample2size[sample])
        rows.append(row)
    return rows, samples

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

def sample_top_clones(clones, minsize, attr='freq'):
    topclone2size = {}
    for clone in clones:
        id = clone.get_vseqj()
        size = clone[attr]
        if size >= minsize:
            if id in topclone2size:
                topclone2size[id] += size
            else:
                topclone2size[id] = size
    return topclone2size

def top_clones(sampledir, minsize, attr='freq'):
    topclone2name2size = {}
    for sam in os.listdir(sampledir):
        samdir = os.path.join(sampledir, sam)
        for vj in os.listdir(samdir):
            if vj == sam:
                continue
            vjfile = os.path.join(samdir, vj)
            clones = pickle.load(gzip.open(vjfile, 'rb'))
            topclone2size = sample_top_clones(clones, minsize, attr)
            for tc, size in topclone2size.iteritems():
                if tc not in topclone2name2size:
                    topclone2name2size[tc] = {sam: size}
                else:
                    topclone2name2size[tc][sam] = size
    return topclone2name2size

def tab_track_clones(clone2rows, groups, outfile):
    f = open(outfile, 'w')
    f.write("#Clone\tIndex\t%s\n" % "\t".join(groups))
    for clone, (rows, indices) in clone2rows.iteritems():
        for i, row in enumerate(rows):
            index = indices[i]
            #rowstr = "\t".join(["%.2e" % c for c in row])
            rowstr = "\t".join(["%.4f" % c for c in row])
            f.write("%s\t%d\t%s\n" % (clone, index, rowstr))
    f.close()

def tab_track_clones_no_matched(clone2rows, groups, outfile):
    f = open(outfile, 'w')
    f.write("#Clone\t%s\n" % "\t".join(groups))
    for clone, (rows, sams) in clone2rows.iteritems():
        f.write("%s" % clone)
        for row in rows:
            f.write("\t%s" % (",".join([str(r) for r in row])))
        f.write("\n")
    f.close()

#====== JobTree Targets ======
class TrackClonesAgg(StatAnalyses):
    '''
    '''
    def __init__(self, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)
    
    def run(self):
        self.load_indir()
        n2o = self.name2obj
        outfile = os.path.join(self.outdir, "trackclones.txt")
        if self.opts.matched:
            tab_track_clones(n2o, self.opts.groups, outfile)
        else:
            tab_track_clones_no_matched(n2o, self.opts.groups, outfile)
        if self.opts.makeplots:
            plotdir = os.path.join(self.outdir, "plots")
            system('mkdir -p %s' % plotdir)
            for clone, (rows, indices) in n2o.iteritems():
                outbase = os.path.join(plotdir, clone)
                if self.opts.matched:
                    tcplot.draw_track_clone_hack(clone, rows, self.opts.groups,
                                        outbase, self.opts)
                    #tcplot.draw_track_clone(clone, rows, self.opts.groups,
                    #                    outbase, self.opts)
                else:
                    tcplot.draw_track_clone_no_matched(clone, rows,
                                         self.opts.groups, outbase, self.opts)

class TrackClone(Target):
    '''track a specific clone
    '''
    def __init__(self, clone, indir, outfile, opts, s2g):
        Target.__init__(self)
        self.clone = clone
        self.indir = indir
        self.outfile = outfile
        self.opts = opts
        self.s2g = s2g

    def run(self):
        items = self.clone.split('_')
        #if len(items) != 3:
        #    print self.clone
        #    print len(items)
        #hack:
        if len(items) > 3:
            items = [items[0], items[-2], items[-1]]
        assert len(items) == 3
        #vj = "%s_%s" % (items[0], items[2])
        v = items[0]
        
        #name2size = clone_get_samples(self.clone, self.indir, vj, 'freq')
        name2size = clone_get_samples(self.clone, self.indir, v, 'freq')
        if self.opts.matched:
            rows, indices = track_clone(self.clone, name2size,
                                        self.opts.groups,
                                        self.opts.group2samples, self.s2g)
            pickle.dump((rows, indices), gzip.open(self.outfile, "wb"))
        else:
            rows, samples = track_clone_no_matched(self.clone, name2size,
                                     self.opts.groups, self.opts.group2samples)
            pickle.dump((rows, samples), gzip.open(self.outfile, "wb"))

class TrackClones(Analysis):
    '''Track a number of input clones
    '''
    def __init__(self, indir, outdir, opts, clones):
        Analysis.__init__(self, indir, outdir, opts)
        self.clones = clones
    
    def run(self):
        opts = self.opts
        sample2group = libcommon.get_val2key_1to1(opts.group2samples)
        global_dir = self.getGlobalTempDir()
        tc_dir = os.path.join(global_dir, "trackclone_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % tc_dir)
        for clone in self.clones:
            outfile = os.path.join(tc_dir, "%s.pickle" % clone)
            self.addChildTarget(TrackClone(clone, self.indir, outfile, opts,
                                           sample2group))
        self.setFollowOnTarget(TrackClonesAgg(tc_dir, self.outdir, opts))

class TrackTopClones(TrackClones):
    '''Set up children jobs to track top clones
    '''
    def __init__(self, indir, outdir, opts):
        topclones = top_clones(indir, opts.track_minfreq)
        TrackClones.__init__(self, indir, outdir, opts, topclones)

