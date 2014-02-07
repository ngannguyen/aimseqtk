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
    return portion, list(common_names)

def group_major_clones(clone2name2size, names, minfreq):
    total = len(names)
    assert total > 0
    majorclones = {}  
    group_names = Set(names)
    for clone, name2size in clone2name2size.iteritems():
        clone_names = Set(name2size.keys())
        common_names = group_names.intersection(clone_names)
        diff_names = clone_names.difference(group_names)
        numsam = len(common_names)
        portion = float(numsam) / total
        if portion >= minfreq:
            #majorclones[clone] = portion
            majorclones[clone] = (numsam, total - numsam, list(common_names))
    return majorclones

def split_c2n2s_by_j(c2n2s):
    j2c2n2s = {}
    for c, n2s in c2n2s.iteritems():
        j = c.rstrip('\n').split('_')[2]
        if j not in j2c2n2s:
            j2c2n2s[j] = {c: n2s}
        else:
            j2c2n2s[j][c] = n2s
    return j2c2n2s

def split_clones_by_j(clones):
    j2clones = {}
    for c, v in clones.iteritems():
        j = c.split('_')[2]
        if j not in j2clones:
            j2clones[j] = {c: v}
        else:
            j2clones[j][c] = v
    return j2clones

class OverlapPairGroupsVj2(Target):
    def __init__(self, mjclones, c2n2s, outfile, names2, opts):
        Target.__init__(self)
        self.mjclones = mjclones
        self.c2n2s = c2n2s
        self.outfile = outfile
        self.names2 = names2
        self.opts = opts

    def run(self):
        f = open(self.outfile, 'w')
        for clone, (present1, absent1, inames) in self.mjclones.iteritems():
            portion1 = float(present1) / (present1 + absent1)
            portion2, onames = clone_group_freq(self.c2n2s, self.names2, clone)
            if portion2 <= self.opts.outgroup:
                present2 = round(portion2*len(self.names2))
                absent2 = len(self.names2) - present2
                tab = [[present1, present2], [absent1, absent2]]
                oddratio, pval = fisher_exact(tab, alternative='greater')
                if pval <= self.opts.pval:
                    p1 = libcommon.pretty_float(portion1)
                    p2 = libcommon.pretty_float(portion2)
                    odd = libcommon.pretty_float(oddratio)
                    p = libcommon.pretty_float(pval)
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (clone, p1, p2,
                                 odd, p, ",".join(inames), ",".join(onames)))
                    #f.write("%s\t%.2e\t%.2e\t%.2e\t%.2e\n" % (clone, present1,
                    #                                 present2, oddratio, pval))
        f.close()
        # remove outfile if there is no significant results
        if os.stat(self.outfile).st_size <= 0:
            system("rm -f %s" % self.outfile)

class OverlapPairGroupsAgg(Target):
    def __init__(self, indir, outfile, header=True):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
        self.header = header

    def run(self):
        if self.header:
            f = open(self.outfile, 'w')
            f.write(("#Clone\tFreq1\tFreq2\tOdd_ratio\tp_value\t" +
                    "Ingroup_samples\tOutgroup_samples\n"))
            f.close()
            
        for file in os.listdir(self.indir):
            filepath = os.path.join(self.indir, file)
            system("cat %s >> %s" % (filepath, self.outfile))
        self.setFollowOnTarget(libcommon.CleanupDir(self.indir))

class OverlapPairGroupsVj(Target):
    def __init__(self, vjfile, c2n2s_file, outfile, names2, opts):
        Target.__init__(self)
        self.vjfile = vjfile
        self.c2n2s_file = c2n2s_file
        self.outfile = outfile
        self.names2 = names2
        self.opts = opts

    def run(self):
        mjclones = pickle.load(gzip.open(self.vjfile, 'rb'))
        c2n2s = pickle.load(gzip.open(self.c2n2s_file, 'rb'))

        # split by j
        j2mjclones = split_clones_by_j(mjclones)
        j2c2n2s = split_c2n2s_by_j(c2n2s)

        tempdir = "%s-jsplit" % self.outfile
        system("mkdir -p %s" % tempdir)
        for j, j_mjclones in j2mjclones.iteritems():
            j_outfile = os.path.join(tempdir, j)
            j_c2n2s = j2c2n2s[j]
            self.addChildTarget(OverlapPairGroupsVj2(j_mjclones, j_c2n2s,
                                            j_outfile, self.names2, self.opts))
        self.setFollowOnTarget(OverlapPairGroupsAgg(tempdir, self.outfile,
                                                    header=False))

class OverlapPairGroups(Target):
    '''Find group_dominant clones of group1
       For each clone, perform fisher exact test
       If matched, track its frequencies over different groups
    '''
    def __init__(self, g1, g2, mjclones_dir, names2, c2n2s_dir, outdir, opts):
        Target.__init__(self)
        self.g1 = g1
        self.g2 = g2
        self.mjclones_dir = mjclones_dir
        self.names2 = names2
        self.c2n2s_dir = c2n2s_dir
        self.outdir = outdir
        self.opts = opts

    def run(self):
        pair = "%s_%s" % (self.g1, self.g2)
        workdir = os.path.join(self.outdir, pair)
        system("mkdir -p %s" % workdir)
        for vj in os.listdir(self.mjclones_dir):
            vj_file = os.path.join(self.mjclones_dir, vj)
            c2n2s_file = os.path.join(self.c2n2s_dir, vj)
            vj_out = os.path.join(workdir, vj)
            self.addChildTarget(OverlapPairGroupsVj(vj_file, c2n2s_file,
                                               vj_out, self.names2, self.opts))
        outfile = os.path.join(self.outdir, "%s.txt" % pair)
        self.setFollowOnTarget(OverlapPairGroupsAgg(workdir, outfile, header=True))
        
class OverlapAllPairGroups(Analysis):
    '''Set up pairwise group analyses
    '''
    def __init__(self, major_clones_dir, outdir, opts, c2n2s_dir):
        Analysis.__init__(self, major_clones_dir, outdir, opts)
        self.c2n2s_dir = c2n2s_dir

    def run(self):
        g2n = self.opts.group2samples
        if g2n and len(g2n) >= 2:
            groups = g2n.keys()
            for g1 in groups:
                #major_clones = self.name2obj[g1]
                g1_major_clones_dir = os.path.join(self.indir, g1) 
                for g2 in groups:
                    if g2 == g1:
                        continue
                    names2 = g2n[g2]
                    self.addChildTarget(OverlapPairGroups(g1, g2,
                                                  g1_major_clones_dir,
                                                  names2, self.c2n2s_dir,
                                                  self.outdir, self.opts))

class GroupClonesVj(Target):
    '''
    '''
    def __init__(self, names, outfile, infile, ingroup):
        Target.__init__(self)
        self.names = names
        self.outfile = outfile
        self.infile = infile
        self.ingroup = ingroup

    def run(self):
        c2n2s = pickle.load(gzip.open(self.infile, 'rb'))
        majorclones = group_major_clones(c2n2s, self.names, self.ingroup)
        if majorclones:
            pickle.dump(majorclones, gzip.open(self.outfile, "wb"))

class GroupClones(Target):
    '''Caculate clones that are present in >= ingroup portion of total
    group samples (major_clones) and clones that are present in <=
    outgroup portion of total group samples (minor_clones)
    '''
    #def __init__(self, names, outfile, clone2name2size, ingroup):
    def __init__(self, names, outdir, c2n2s_dir, ingroup):
        Target.__init__(self)
        self.names = names
        self.outdir = outdir
        #self.clone2name2size = clone2name2size
        self.c2n2s_dir = c2n2s_dir
        self.ingroup = ingroup

    def run(self):
        for vj in os.listdir(self.c2n2s_dir):
            c2n_file = os.path.join(self.c2n2s_dir, vj)
            outfile = os.path.join(self.outdir, vj)
            self.addChildTarget(GroupClonesVj(self.names, outfile,
                                              c2n_file, self.ingroup))
        #self.setFollowOnTarget(GroupClonesAgg(workdir, self.outfile))

class Overlap2(Target):
    def __init__(self, c2n2s_dir, outdir, opts):
        Target.__init__(self)
        self.c2n2s_dir = c2n2s_dir
        self.outdir = outdir
        self.opts = opts

    def run(self):
        #clone2sample2size = {}
        #objs = libcommon.load_pickledir(self.c2n2s_dir)
        #for c2n2s in objs:
        #    clone2sample2size.update(c2n2s)

        g2n = self.opts.group2samples
        if g2n and len(g2n) >= 2:
            global_dir = self.getGlobalTempDir()
            o_dir = os.path.join(global_dir, "overlap_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
            system("mkdir -p %s" % o_dir)
            for group, names in g2n.iteritems():  # get major & minor clones
                #groupfile = os.path.join(o_dir, "%s.pickle" % group)  # pickle file
                groupdir = os.path.join(o_dir, group)
                system('mkdir -p %s' % groupdir)
                self.addChildTarget(GroupClones(names, groupdir,
                                         self.c2n2s_dir, self.opts.ingroup))
            self.setFollowOnTarget(OverlapAllPairGroups(o_dir,
                               self.outdir, self.opts, self.c2n2s_dir))

class Overlap(Analysis):
    '''Set up children jobs to compute clones that are dominantly
    present in a specific group
    '''
    def __init__(self, indir, outdir, opts):
        Analysis.__init__(self, indir, outdir, opts)

    def run(self):
        self.logToMaster("Overlap\n")
        global_dir = self.getGlobalTempDir()
        workdir = os.path.join(global_dir, "overlap_temp")
        system('mkdir -p %s' % workdir)
        sizetype = 'freq'
        self.addChildTarget(statcommon.GetClone2Samples(self.indir, workdir,
                                                        sizetype))
        self.setFollowOnTarget(Overlap2(workdir, self.outdir, self.opts))



