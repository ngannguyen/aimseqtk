#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Get clones that are present in >= minNumInGroup and in <= maxNumOutGroup
Perform Fisher exact test for each of those clones
Bon ferroni correction
'''

import os
import sys
import gzip
import cPickle as pickle
from sets import Set
from optparse import OptionGroup

from scipy.stats import fisher_exact
from numpy import median
from jobTree.scriptTree.target import Target
from sonLib.bioio import system

from aimseqtk.lib.statcommon import SampleStat
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
from aimseqtk.lib.common import StatAnalyses
from aimseqtk.lib.common import Analysis
import aimseqtk.lib.statcommon as statcommon


def add_overlap_options(parser):
    group = OptionGroup(parser, "Group_dominant_clones options")
    group.add_option('--min_in_group', dest='ingroup', default=1.0,
                     type='float',
                     help=('Minimum proportion of in_group samples. ' +
                           'Default=%default. Ranges from 0.0 to 1.0'))
    group.add_option('--max_out_group', dest='outgroup', default=0.0,
                     type='float',
                     help=('Max proportion of out_group samples. ' +
                           'Default=%default. Ranges from 0.0 to 1.0'))
    parser.add_option_group(group)
                    
def clone_group_freq(clone2name2size, names, clone):
    # proprotion of group's sample that have "clone"
    assert len(names) > 0
    if clone not in clone2name2size:
        return 0.0
    all_names = Set(clone2name2size[clone].keys())
    group_names = Set(names)
    common_names = group_names.intersection(all_names)
    portion = float(len(common_names)) / len(names)
    return portion

def group_major_clones(clone2name2size, names, minfreq):
    total = len(names)
    assert total > 0
    majorclones = {}  
    group_names = Set(names)
    for clone, name2size in clone2name2size.iteritems():
        clone_names = Set(name2size.keys())
        common_names = group_names.intersection(clone_names)
        portion = float(len(common_names)) / total
        if portion >= minfreq:
            majorclones[clone] = portion
    return majorclones

class OverlapPairGroups(Target):
    '''Find group_dominant clones of group1
       For each clone, perform fisher exact test
       If matched, track its frequencies over different groups
    '''
    def __init__(self, g1, g2, major_clones, names2, clone2name, outdir, opts):
        Target.__init__(self)
        self.g1 = g1
        self.g2 = g2
        self.major_clones = major_clones
        self.names2 = names2
        self.clone2name2size = clone2name
        self.outdir = outdir
        self.opts = opts

    def run(self):
        pair = "%s_%s" % (self.g1, self.g2)
        signi_clones = []
        outfile = os.path.join(self.outdir, "%s.txt" % pair)
        f = open(outfile, 'w')
        f.write("#Clone\tFreq1\tFreq2\tOdd_ratio\tp_value\n")
        
        for clone, portion in self.major_clones.iteritems():
            present2 = clone_group_freq(self.clone2name2size, self.names2, clone)
            if present2 <= self.opts.outgroup:
                present1 = self.major_clones[clone]
                absent1 = 1.0 - present1
                absent2 = 1.0 - present2
                tab = [[present1, present2], [absent1, absent2]] 
                oddratio, pval = fisher_exact(tab, alternative='greater')
                if pval <= self.opts.pval:
                    f.write("%s\t%.2e\t%.2e\t%.2e\t%.2e\n" % (clone, present1,
                                                     present2, oddratio, pval))
                signi_clones.append(clone)
                #if self.opts.matched:
                #    draw_clone_freq_over_diff_group
        f.close()
        # remove outfile if there is no significant results
        if os.stat(outfile).st_size <= 0:
            system("rm -f %s" % outfile)
             
class OverlapAllPairGroups(StatAnalyses):
    '''Set up pairwise group analyses
    '''
    def __init__(self, indir, outdir, opts, c2n2s):
        StatAnalyses.__init__(self, indir, outdir, opts)
        self.clone2name2size = c2n2s

    def run(self):
        self.load_indir()
        g2n = self.opts.group2samples
        if g2n and len(g2n) >= 2:
            groups = g2n.keys()
            for g1 in groups:
                names1 = g2n[g1]
                major_clones = self.name2obj[g1]
                for g2 in groups:
                    if g2 == g1:
                        continue
                    names2 = g2n[g2]
                    self.addChildTarget(OverlapPairGroups(g1, g2, major_clones,
                                                  names2, self.clone2name2size,
                                                  self.outdir, self.opts))

class GroupClones(Target):
    '''Caculate clones that are present in >= ingroup portion of total
    group samples (major_clones) and clones that are present in <=
    outgroup portion of total group samples (minor_clones)
    '''
    def __init__(self, names, outfile, clone2name2size, ingroup):
        Target.__init__(self)
        self.names = names
        self.outfile = outfile
        self.clone2name2size = clone2name2size
        self.ingroup = ingroup

    def run(self):
        majorclones = group_major_clones(self.clone2name2size, self.names,
                                                                  self.ingroup)
        pickle.dump(majorclones, gzip.open(self.outfile, "wb"))

class Overlap(Analysis):
    '''Set up children jobs to compute clones that are dominantly
    present in a specific group
    '''
    def __init__(self, samples, outdir, opts):
        Analysis.__init__(self, samples, outdir, opts)

    def run(self):
        sizetype = 'freq'
        clone2sample2size = statcommon.get_clone2samples(self.samples, sizetype)
        g2n = self.opts.group2samples
        if g2n and len(g2n) >= 2:
            global_dir = self.getGlobalTempDir()
            o_dir = os.path.join(global_dir, "overlap_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
            system("mkdir -p %s" % o_dir)
            for group, names in g2n.iteritems():  # get major & minor clones
                groupfile = os.path.join(o_dir, "%s.pickle" % group)  # pickle file
                self.addChildTarget(GroupClones(names, groupfile,
                                         clone2sample2size, self.opts.ingroup))
            self.setFollowOnTarget(OverlapAllPairGroups(o_dir,
                               self.outdir, self.opts, clone2sample2size))



