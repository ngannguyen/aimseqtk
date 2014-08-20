#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Compute similarity between pair of samples
1/ Similarity indices include: 1/number of overlap clones 2/ chao 3/ horn, etc
For each index, create a heatmap with tree (rows=samples, cols=samples, cells=index)
2/ Test for similarity differences: for each pair of groups, compare: within_group 
samples and between_group samples
3/ Plot: # shared clones vs # samples, one line/group
4/ Plot: # shared clones vs n1*n2 (for each samplee; one line/other_sample)
5/ Test for the presence/ absence of clone in different groups

Note: R vegan package support similarity indices:
"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower",
"morisita", "horn", "mountford", "raup" , "binomial" or "chao"
Will use: ['horn', 'chao']
'''

import os
import sys
import numpy as np
import gzip
#import marshal as pickle
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
from aimseqtk.lib.common import Analysis
from aimseqtk.lib.common import StatAnalyses
import aimseqtk.lib.tabcommon as tabcommon
import aimseqtk.lib.statcommon as statcommon
import aimseqtk.lib.drawcommon as drawcommon
import aimseqtk.src.properties.similarity_plot as simiplot


class SimilarityStat(statcommon.PairStat):
    def __init__(self):
        statcommon.PairStat.__init__(self)
        self.numshare = None
        self.chao = None
        self.horn = None

def extract_c2n2s(names, c2n2s_dir):
    c2n2s_pair = {}  # clones that are present in >=1 sample in "names"
    for v in os.listdir(c2n2s_dir):
        c2n_file = os.path.join(c2n2s_dir, v)
        c2n2s = pickle.load(gzip.open(c2n_file, 'rb'))
        for c, n2s in c2n2s.iteritems():
            n2s_pair = {n: n2s[n] for n in names if n in n2s}
            if n2s_pair:
                c2n2s_pair[c] = n2s_pair
    return c2n2s_pair
    
def pair_similarity(sample1, sample2, c2n2s, attrs):
    total = len(c2n2s)
    # get count vectors:
    vec1 = []
    vec2 = []
    numshare = 0
    for clone, sample2size in c2n2s.iteritems():
        v1 = 0.0
        if sample1.name in sample2size:
            v1 = sample2size[sample1.name]
        vec1.append(v1)
        v2 = 0.0
        if sample2.name in sample2size:
            v2 = sample2size[sample2.name]
        vec2.append(v2)
        if v1 > 0 and v2 > 0:
            numshare += 1

    # get similarity indices:        
    stat = SimilarityStat()
    stat.set_sample_info(sample1, sample2)
    
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    vegan = importr("vegan")
    vegdist = vegan.vegdist

    vec = vec1 + vec2
    rvec = ro.FloatVector(vec)
    rmatrix = ro.r['matrix'](rvec, nrow=2, byrow=True)
    for attr in attrs:
        if attr == 'numshare':
            stat[attr] = numshare
        else:
            dissimilarity = vegdist(rmatrix, method=attr)
            stat[attr] = 1.0 - dissimilarity[0]
    return stat

def check_similarity_options(options):
    indices = []
    valid_indices = ['numshare', 'chao', 'horn']
    if options.similarity:
        items = options.similarity.split(',')
        for item in items:
            if item in valid_indices:
                indices.append(item)
    options.similarity = indices            
    
def add_similarity_options(parser):
    group = OptionGroup(parser, "Similarity options")
    group.add_option('--similarity', dest='similarity', default='numshare',
                     help=('Comma separated list of similarity indices. ' +
                           'Default=%default. Valid options are: [numshare' +
                           ',chao,horn.'))
    parser.add_option_group(group)

class SimilarityHeatmap(Target):
    ''' Make heatmap table and heatmap plot for a specific similarity index
    '''
    def __init__(self, samplenames, pair2stat, outbase, attr, opts):
        Target.__init__(self)
        self.names = samplenames
        self.pair2stat = pair2stat
        self.outbase = outbase
        self.attr = attr
        self.opts = opts

    def run(self):
        # Get the matrix
        rows = statcommon.pair_matrix(self.names, self.names, self.pair2stat,
                                                                     self.attr)
        # Write to table
        tabfile = "%s.txt" % self.outbase
        tabcommon.matrix_table(self.names, rows, tabfile)
        # Make heatmap
        if self.opts.makeplots:
            plotfile = "%s.pdf" % self.outbase
            drawcommon.draw_heatmap(self.names, self.names, rows, plotfile)

def islarge(vecs):
    minsam = 20
    large = True
    for v in vecs:
        if len(v) < 20:
            large = False
            break
    return large

def table_group_pairwise_similarity(g1, g2, vec11, vec12, vec22, outfile, testtype=None):
    mean_11, std_11 = statcommon.vec_mean_std(vec11)
    mean_12, std_12 = statcommon.vec_mean_std(vec12)
    mean_22, std_22 = statcommon.vec_mean_std(vec22)
    
    f = open(outfile, 'w')
    f.write("Categories\tpval\ttval\tmean1 +/- std1\tmean2 +/- std2\n")
    # Compare 11 and 12
    if testtype is None:
        testtype = statcommon.get_test_type2(islarge([vec11, vec12]), False)
    f.write("%s, %s_%s\t" % (g1, g1, g2))
    t_11_12, p_11_12 = statcommon.ttest_pair(vec11, vec12, testtype=testtype) 
    f.write("%.3f\t%.2e\t%.3f +/- %.3f\t%.3f +/- %.3f\n" % (t_11_12, p_11_12,
                                            mean_11, std_11, mean_12, std_12))
    # Compare 12 and 22
    if testtype is None:
        testtype = statcommon.get_test_type2(islarge([vec12, vec22]), False)
    f.write("%s_%s, %s\t" % (g1, g2, g2))
    t_12_22, p_12_22 = statcommon.ttest_pair(vec12, vec22, testtype=testtype) 
    f.write("%.3f\t%.2e\t%.3f +/- %.3f\t%.3f +/- %.3f\n" % (t_12_22, p_12_22,
                                            mean_12, std_12, mean_22, std_22))
    # Compare 11 and 22
    if testtype is None:
        testtype = statcommon.get_test_type2(islarge([vec11, vec22]), False)
    f.write("%s, %s\t" % (g1, g2))
    t_11_22, p_11_22 = statcommon.ttest_pair(vec11, vec22, testtype=testtype) 
    f.write("%.3f\t%.2e\t%.3f +/- %.3f\t%.3f +/- %.3f\n" % (t_11_22, p_11_22,
                                            mean_11, std_11, mean_22, std_22))
    f.close()

class SimilarityPairGroups(Target):
    '''Within and between group comparisons for a specific similarity index
    '''
    def __init__(self, pair2stat, g1, names1, g2, names2, attr, out, opts):
        Target.__init__(self)
        self.pair2stat = pair2stat
        self.g1 = g1
        self.names1 = names1
        self.g2 = g2
        self.names2 = names2
        self.attr = attr
        self.outbase = out
        self.opts = opts

    def run(self):
        vec11 = statcommon.pair_vec(self.names1, self.names1,
                                        self.pair2stat, self.attr)
        vec22 = statcommon.pair_vec(self.names2, self.names2,
                                        self.pair2stat, self.attr)
        vec12 = statcommon.pair_vec(self.names1, self.names2,
                                        self.pair2stat, self.attr)
        # draw boxplots
        plotfmt = self.opts.plotformat
        plotfile = "%s" % self.outbase
        simiplot.draw_group_pairwise_similarity(self.g1, self.g2, vec11,
                           vec12, vec22, plotfile, plotfmt, self.opts.dpi)
        # ttests 
        tabfile = "%s.txt" % self.outbase
        table_group_pairwise_similarity(self.g1, self.g2, vec11, vec12, vec22,
                                        tabfile)

class SimilarityAnalyses(StatAnalyses):
    ''' Set children jobs to make heatmaps for each similarity index
        and Compare each pair of groups
    '''
    def __init__(self, samplenames, indir, outdir, opts):
        StatAnalyses.__init__(self, indir, outdir, opts)
        self.names = samplenames

    def run(self):
        self.logToMaster("SimilarityAnalyses\n")
        self.load_indir()
        opts = self.opts
        pair2stat = self.name2obj
        outdir = os.path.join(self.outdir, "similarity")
        system("mkdir -p %s" % outdir)

        # Heatmap plots and tables of each similarity index
        heatdir = os.path.join(outdir, "heatmap")
        system("mkdir -p %s" % heatdir)
        for attr in opts.similarity:
            outbase = os.path.join(heatdir, attr)
            self.addChildTarget(SimilarityHeatmap(self.names, pair2stat,
                                                  outbase, attr, opts))
        # Test similarity differences for each pair of groups
        g2n = opts.group2samples
        if g2n:
            cmpdir = os.path.join(outdir, "group_comparisons")
            system("mkdir -p %s" % cmpdir)
            groups = g2n.keys()
            for i in xrange(0, len(groups) - 1):
                g1 = groups[i]
                names1 = g2n[g1]
                for j in xrange(i + 1, len(groups)):
                    g2 = groups[j]
                    names2 = g2n[g2]
                    for attr in opts.similarity:
                        out = os.path.join(cmpdir, "%s_%s_%s" % (g1, g2, attr))
                        self.addChildTarget(SimilarityPairGroups(pair2stat,
                                      g1, names1, g2, names2, attr, out, opts))
        self.setFollowOnTarget(libcommon.CleanupDir(self.indir))

class PairSimilarity(Target):
    '''Compute similarity indices for a specific pair of sample
    Pickle the outputs to outfile
    '''
    def __init__(self, sample1, sample2, c2n2s_dir, outfile, opts):
        Target.__init__(self)
        self.sample1 = sample1
        self.sample2 = sample2
        self.c2n2s_dir = c2n2s_dir
        self.outfile = outfile
        self.opts = opts

    def run(self):
        self.logToMaster("PairSimilarity\n")
        n1 = self.sample1.name
        n2 = self.sample2.name
        c2n2s_pair = extract_c2n2s([n1, n2], self.c2n2s_dir)
        stat = pair_similarity(self.sample1, self.sample2, c2n2s_pair,
                               self.opts.similarity)
        pickle.dump(stat, gzip.open(self.outfile, 'wb'))

class Similarity2(Target):
    def __init__(self, c2n2s_dir, sams, indir, outdir, opts):
        Target.__init__(self)
        self.c2n2s_dir = c2n2s_dir
        self.sams = sams
        self.outdir = outdir
        self.indir = indir
        self.opts = opts

    def run(self):
        global_dir = self.getGlobalTempDir()
        s_dir = os.path.join(global_dir, "similarity_%s" %
                                     os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % s_dir)
        numsam = len(self.sams)
        for i1 in xrange(numsam - 1):
            sam1 = self.sams[i1]
            sam1_file = os.path.join(self.indir, sam1, sam1)
            sample1 = pickle.load(gzip.open(sam1_file, 'rb'))
            for i2 in xrange(i1 + 1, numsam):
                sam2 = self.sams[i2]
                sam2_file = os.path.join(self.indir, sam2, sam2)
                sample2 = pickle.load(gzip.open(sam2_file, 'rb'))
                outfile = os.path.join(s_dir, "%s_%s.pickle" % (sam1, sam2))
                self.addChildTarget(PairSimilarity(sample1, sample2,
                                                   self.c2n2s_dir,
                                                   outfile, self.opts))
        self.setFollowOnTarget(SimilarityAnalyses(self.sams, s_dir,
                                                  self.outdir, self.opts))
        
class Similarity(Analysis):
    ''' Set up children jobs to compute similarity indices for each pair
    of samples
    '''
    def __init__(self, indir, outdir, opts):
        Analysis.__init__(self, indir, outdir, opts)

    def run(self):
        self.logToMaster("Similarity\n")
        sams = os.listdir(self.indir)
        numsam = len(sams)
        if numsam < 2:
            return
        workdir = "similarity_temp"
        system("mkdir -p %s" % workdir)
        sizetype = 'count'
        if self.opts.normalize:
            sizetype = 'normcount'
        self.addChildTarget(statcommon.GetClone2Samples(self.indir, workdir,
                                                        sizetype))
        self.setFollowOnTarget(Similarity2(workdir, sams, self.indir,
                                           self.outdir, self.opts))


