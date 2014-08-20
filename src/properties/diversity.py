#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Perform rarefaction analyses and compute diversity indices
Plot: xaxis: # sampling size; yaxis: diversity index. One curve/sample
'''

import os
import sys
import time
import numpy as np
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system
from sonLib.bioio import logger

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
from aimseqtk.lib.common import Analysis
from aimseqtk.lib.common import StatAnalyses
import aimseqtk.lib.tabcommon as tabcommon
import aimseqtk.lib.statcommon as statcommon
from aimseqtk.lib.statcommon import SampleStat
import aimseqtk.src.properties.diversity_plot as dvplot
import aimseqtk.src.properties.rarefaction_plot as rfplot


class SampleDiversityStat(SampleStat):
    '''Different indices commpute from samplings of a specific size of
       a sample
    '''
    def __init__(self):
        SampleStat.__init__(self)
        #self.numclone = None
        # R vegan: diversity
        self.simpson = None 
        self.invsimpson = None
        self.shannon = None
        self.fisher_alpha = None

        # Standard deviation
        self.numclone_std = None
        self.simpson_std = None 
        self.invsimpson_std = None
        self.shannon_std = None
        self.fisher_alpha_std = None

def rf_diversity_table(name2size2sampling, outfile, index):
    # Rows = Sampling_Sizes, Columns = Samples
    names = sorted(name2size2sampling.keys())
    sizes = []
    for name, size2sampling in name2size2sampling.iteritems():
        for size in size2sampling:
            if size not in sizes:
                sizes.append(size)
    sizes = sorted(sizes)

    f = open(outfile, 'w')
    f.write("#Sample_size\t%s\n" % "\t".join(names))
    for size in sizes:
        f.write("%d" % size)
        for name in names:
            if size not in name2size2sampling[name]:
                f.write("\tNA")
            else:
                sampling = name2size2sampling[name][size]
                if index not in sampling.getitems():
                    f.write("\tNA")
                else:
                    f.write("\t%.3f" % sampling[index])
        f.write("\n")
    f.close()

def sample_get_clone_attrs(name, samdir, attr):
    # return a list of attr of all clones of sample
    vec = []
    for vj in os.listdir(samdir):
        if vj == name:
            continue
        vjfile = os.path.join(samdir, vj)
        clones = pickle.load(gzip.open(vjfile, "rb"))
        vjvec = [c[attr] for c in clones]
        vec.extend(vjvec)
    return vec

def sample_sampling_diversity(sample, samdir, args=None):
    # Sampling sample to "size, and compute the "index" of the subsample
    assert len(args) == 3
    size = args[0]
    indices = args[1]
    workdir = args[2]
    if size is None:  # no sampling, compute diversity of original sample
        subsampledir = samdir
    else:
        subsampledir = os.path.join(workdir, "subsample")
        system("mkdir -p %s" % subsampledir)
        libsample.sampling(sample, samdir, subsampledir, (size,))
    
    sampling = SampleDiversityStat()
    sampling.set_sample_info(sample)
    counts = sample_get_clone_attrs(sample.name, subsampledir, "count")

    #import rpy2.rinterface as rinterface
    #rinterface.set_initoptions(('rpy2', '--no-save', '--no-restore'))
    #rinterface.initr()

    import rpy2.robjects as robjs
    from rpy2.robjects.packages import importr
    vegan = importr("vegan")
    rcounts = robjs.IntVector(counts)
    for index in indices:
        if index == 'numclone':
            sampling[index] = len(counts)
        elif index == 'fisher_alpha':
            rfisher = vegan.fisher_alpha(rcounts)
            sampling[index] = rfisher[0]
        else:
            rval = vegan.diversity(rcounts, index)
            sampling[index] = rval[0]
    return sampling

def sample_avr_sampling_diversity(samplings, indices):
    assert len(samplings) >= 1
    avrsampling = SampleDiversityStat()
    s1 = samplings[0]
    avrsampling.set_sample_info2(s1.name, s1.group, s1.color, s1.marker)

    stds = ["%s_std" % i for i in indices]
    for i, index in enumerate(indices):
        vals = [sampling[index] for sampling in samplings]
        avrsampling[index] = np.mean(vals)
        avrsampling[stds[i]] = np.std(vals)
    return avrsampling

def sample_rf_sizes(sample_size, bin=None, rf_sizes=None):
    # Return the list of sampling size for the sample
    sizes = []
    if rf_sizes:
        for s in rf_sizes:
            if s <= sample_size:
                sizes.append(s)
    elif bin and bin <= sample_size:
        sizes = range(bin, sample_size + 1, bin)
    
    if len(sizes) == 0:
        sizes = [sample_size]
    return sizes

def diversity_ttest_text(index2obj, outfile, pcutoff):
    # obj: (pair2(tval, pval), group2(mean, std))
    f = open(outfile, 'w')
    f.write("#Index\tSample_Pair\tt_value\tp_value\tMean1\tStd1\tMean2\tStd2\n")
    for index in sorted(index2obj.keys()):
        obj = index2obj[index]
        assert len(obj) == 2
        pair2stat = obj[0]
        group2mean = obj[1]
        for pair in sorted(pair2stat.keys()):
            tval, pval = pair2stat[pair]
            if pval > pcutoff:  # does not pass pcutoff
                continue
            groups = pair.split('_')
            assert len(groups) == 2
            (mean1, std1) = group2mean[groups[0]]
            (mean2, std2) = group2mean[groups[1]]
            vals = [tval, pval, mean1, std1, mean2, std2]
            valstrs = []
            for v in vals:
                valstrs.append(libcommon.pretty_float(v))
            f.write("%s\t%s\t%s\n" % (index, pair, '\t'.join(valstrs)))
    f.close()

def check_rarefaction_options(parser, options):
    if options.rf_sizes:
        options.rf_sizes = [long(s) for s in options.rf_sizes.split(",")]
    indices = []
    valid_indices = ['numclone', 'simpson', 'invsimpson', 'shannon',
                     'fisher_alpha']
    if options.diversity:
        items = options.diversity.split(',')
        for item in items:
            if item not in valid_indices:
                raise ValueError("Unknown diversity index: %s" % item)
            indices.append(item)
    options.diversity = indices

def add_rarefaction_options(parser):
    group = OptionGroup(parser, "Rarefaction analyses options")
    group.add_option('--rf_num_sampling', type='int', default=1,
                     help=('Number of samplings performed for each sampling' +
                           'size. Default=%default'))
    group.add_option('--bin', dest='bin', type='int', default=50000,
                     help=('Increment of sampling size for rarefaction ' +
                           'analyses. Default=%default, i.e sampling sizes ' +
                           'are 50k, 100k,150k, etc.'))
    group.add_option('--rf_sizes', dest='rf_sizes', 
                     help=('Optional. Comma separated string of sampling ' +
                        'sizes. If specified, will ignore options "bin".'))
    group.add_option('--diversity', dest='diversity',
                     default='numclone,shannon,fisher_alpha',
                     help=('Comma separated list of diversity indices. ' +
                           'Default=%default. Valid options are: [numclone' +
                           ',simpson,invsimpson,shannon,fisher_alpha].'))
    parser.add_option_group(group)

#============ PARALLELIZE =============
class AvrSamplingDiversity(Target):
    def __init__(self, indir, indices, outdir, workdir):
        Target.__init__(self)
        self.indir = indir
        self.indices = indices
        self.outdir = outdir
        self.workdir = workdir

    def run(self):
        system("rm -Rf %s" % self.workdir)  # cleanup
        for sizedir in os.listdir(self.indir):
            sizedirpath = os.path.join(self.indir, sizedir)
            samplings = libcommon.load_pickledir(sizedirpath)
            avrsampling = sample_avr_sampling_diversity(samplings, 
                                                        self.indices)
            outfile = os.path.join(self.outdir, "%s.pickle" % sizedir)
            pickle.dump(avrsampling, gzip.open(outfile, "wb"))
        self.setFollowOnTarget(libcommon.CleanupDir(self.indir))

class SampleDiversityRarefaction(Target):
    '''Perform rarefaction analyses on the specific sample
    "numsampling" of times
    '''
    def __init__(self, outdir, sample, samdir, sizes, numsampling, indices,
                                                           avrdir, workdir):
        Target.__init__(self)
        self.outdir = outdir
        self.sample = sample
        self.samdir = samdir
        self.sizes = sizes
        self.numsampling = numsampling
        self.indices = indices
        self.avrdir = avrdir
        self.workdir = workdir

    def run(self):
        for size in self.sizes:
            sizedir = os.path.join(self.outdir, "%d" % size)
            system("mkdir -p %s" % sizedir)
            for i in xrange(self.numsampling):
                outfile = os.path.join(sizedir, "%d" % i)
                workdir = os.path.join(self.workdir, str(size), str(i))
                system("mkdir -p %s" % workdir)
                self.addChildTarget(libsample.SampleAnalysis(self.sample,
                                     self.samdir, outfile,
                                     sample_sampling_diversity,
                                     size, self.indices, workdir))
        self.setFollowOnTarget(AvrSamplingDiversity(self.outdir, self.indices,
                                                    self.avrdir, self.workdir))

class DiversityRarefactionSummary(Target):
    # Print summary tables and draw plots
    def __init__(self, indir, sampledir, options, workdir):
        Target.__init__(self)
        self.indir = indir   # avrdir: avrdir/sample/size.pickle
        self.sampledir = sampledir  # samples/sample/vj
        self.options = options
        self.workdir = workdir

    def run(self):
        system("rm -Rf %s" % self.workdir)
        opts = self.options
        indices = opts.diversity
        name2size2sampling = {}
        size2name2sampling = {}
        for name in os.listdir(self.indir):
            sampledir = os.path.join(self.indir, name)
            name2size2sampling[name] = {}
            for file in os.listdir(sampledir):
                size = long(file.split(".")[0])
                filepath = os.path.join(sampledir, file)
                sampling = pickle.load(gzip.open(filepath, "rb"))
                name2size2sampling[name][size] = sampling

                if size not in size2name2sampling:
                    size2name2sampling[size] = {name: sampling}
                else:
                    size2name2sampling[size][name] = sampling

        outdir = os.path.join(opts.outdir, "diversity")
        txtdir = os.path.join(outdir, "txt")
        system("mkdir -p %s" % txtdir)
        texdir = os.path.join(outdir, "tex")
        system("mkdir -p %s" % texdir)
        pdfdir = os.path.join(outdir, "pdf")
        if opts.makeplots:
            system("mkdir -p %s" % pdfdir)
        
        # Summary table for each index:
        groups = None
        g2s = opts.group2samples
        if g2s:
            groups = opts.group2samples.keys()

        for index in indices:
            rftabfile = os.path.join(txtdir, "rf_%s.txt" % index)
            rf_diversity_table(name2size2sampling, rftabfile, index)
            if opts.makeplots:
                rfplotfile = os.path.join(pdfdir, "rf_%s" % index)
                rfplot.draw_rarefaction(name2size2sampling, groups, index,
                                        rfplotfile, opts.plotformat, opts.dpi)
        # Summary table for each sampling size:
        for size, name2sampling in size2name2sampling.iteritems():
            g2avr = {}
            if g2s:
                g2avr = libcommon.get_group_avr(name2sampling, g2s)
            txtfile = os.path.join(txtdir, "diversity_%d.txt" % size)
            tabcommon.table(name2sampling, txtfile, indices, g2avr, g2s)
            texfile = os.path.join(texdir, "diversity_%d.tex" % size)
            tabcommon.table(name2sampling, texfile, indices, g2avr, g2s, True)
        self.setFollowOnTarget(libcommon.CleanupDir(self.indir))

class DiversityRarefaction(Target):
    # diversity_avr/sample/size.pickle
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        global_dir = self.getGlobalTempDir()
        rf_dir = os.path.join(global_dir, 'rf')
        outdir = os.path.join(rf_dir, "diversity")
        system("mkdir -p %s" % outdir)
        avrdir = os.path.join(rf_dir, "diversity_avr")
        system("mkdir -p %s" % avrdir)
        workdir = os.path.join(rf_dir, "temp")
        system("mkdir -p %s" % workdir)

        bin = self.options.bin
        rf_sizes = self.options.rf_sizes
        for name in os.listdir(self.indir):
            samdir = os.path.join(self.indir, name)
            sam_outdir = os.path.join(outdir, name)
            avr_outdir = os.path.join(avrdir, name)
            sam_workdir = os.path.join(workdir, name)
            system("mkdir -p %s" % sam_outdir)
            system("mkdir -p %s" % avr_outdir)
            system("mkdir -p %s" % sam_workdir)
            sample = pickle.load(gzip.open(os.path.join(samdir, name), "rb"))
            sizes = sample_rf_sizes(sample.size, bin, rf_sizes)
            self.addChildTarget(SampleDiversityRarefaction(sam_outdir, sample,
                             samdir, sizes, self.options.rf_num_sampling,
                             self.options.diversity, avr_outdir, sam_workdir))
        self.setFollowOnTarget(DiversityRarefactionSummary(avrdir,
                                   self.indir, self.options, workdir))

class DiversityTtest(Target):
    #Outfile: index.pickle (pair2tv, group2meanstd)
    def __init__(self, outdir, index, g2n, name2stat, matched, plotfmt,
                                                                      plotdir):
        Target.__init__(self)
        self.outdir = outdir
        self.index = index
        self.g2n = g2n
        self.name2stat = name2stat
        self.matched = matched
        self.plotfmt = plotfmt
        self.plotdir = plotdir

    def run(self):
        testtype = statcommon.get_test_type(self.g2n, self.matched)
        p2stat, g2mean = statcommon.ttest_allpairs(self.g2n, self.name2stat,
                                                   self.matched, self.index,
                                                   testtype=testtype)
        picklefile = os.path.join(self.outdir, "%s.pickle" % self.index)
        pickle.dump((p2stat, g2mean), gzip.open(picklefile, 'wb'))
        if self.plotfmt:
            plotfile = os.path.join(self.plotdir, self.index)
            if not self.matched:
                dvplot.draw_diversity_plot(self.g2n, self.name2stat,
                                           self.index, plotfile, self.plotfmt)
            else:
                dvplot.draw_diversity_plot_hacktimeseries(self.g2n,
                           self.name2stat, self.index, plotfile, self.plotfmt)

class DiversityTtestSummary(StatAnalyses):
    '''Print table summary of the ttest statistics
       Indir: index.pickle (pair2tv, group2meanstd)
       Outtable: cols: Index, PairSamples, tval, pval, mean1 +- std1,
       mean2 +- std2;  rows: each index, each sample pair
    '''
    def __init__(self, indir, outdir, pval=1.0):
        StatAnalyses.__init__(self, indir, outdir)
        self.pval = pval

    def run(self):
        self.load_indir()
        index2obj = self.name2obj
        txtfile = os.path.join(self.outdir, "diversity_ttests.txt")
        diversity_ttest_text(index2obj, txtfile, self.pval)

class DiversitySummary(StatAnalyses):
    #  table of diversity indices
    def __init__(self, sampledir, outdir, indices, group2samples, matched,
                                                plotfmt=None, pval=1.0):
        StatAnalyses.__init__(self, sampledir, outdir)
        self.indices = indices
        self.group2samples = group2samples
        self.matched = matched
        self.plotfmt = plotfmt
        self.pval = pval

    def run(self):
        self.load_indir()
        name2sampling = self.name2obj
        # Print out summary table
        g2s = self.group2samples
        g2avr = {}
        if g2s:
            g2avr = libcommon.get_group_avr(name2sampling, g2s)
        txtfile = os.path.join(self.outdir, "diversity.txt")
        tabcommon.table(name2sampling, txtfile, self.indices, g2avr, g2s)
        texfile = os.path.join(self.outdir, "diversity.tex")
        tabcommon.table(name2sampling, texfile, self.indices, g2avr, g2s, True)

        # For each diversity index, each pair of groups, perform ttest
        # and draw plot
        if g2s and len(g2s) >= 2:
            outdir = os.path.join(self.outdir, "group_comparisons")
            system("mkdir -p %s" % outdir)
            global_dir = os.path.join(self.getGlobalTempDir(), "ttests")
            system("mkdir -p %s" % global_dir)
            for index in self.indices:
                self.addChildTarget(DiversityTtest(global_dir, index, g2s,
                            name2sampling, self.matched, self.plotfmt, outdir))
            self.setFollowOnTarget(DiversityTtestSummary(global_dir, outdir,
                                                         self.pval))

class Diversity(Analysis):
    '''Calculate diversity indices for input samples. No sampling.
       Table: Rows=Samples; Cols=Indices
       If group2samples is provided: perform statistic test comparing
       pair of groups. If plot is True, draw plot: boxplots of: 
       xaxis: groups, yaxis: diversity indices
    '''
    def __init__(self, indir, outdir, indices, group2samples, matched,
                 plotfmt=None, pval=1.0):
        Analysis.__init__(self, indir, outdir)
        self.indices = indices
        self.group2samples = group2samples
        self.matched = matched
        self.plotfmt = plotfmt
        self.pval = pval

    def run(self):
        size = None  # no sampling
        global_dir = self.getGlobalTempDir()

        d_dir = os.path.join(global_dir, "diversity_%s" %
                                    os.path.basename(self.outdir.rstrip('/')))
        system("mkdir -p %s" % d_dir)

        for sam in os.listdir(self.indir):
            samdir = os.path.join(self.indir, sam)
            sample = pickle.load(gzip.open(os.path.join(samdir, sam), 'rb'))
            outfile = os.path.join(d_dir, "%s.pickle" % sam)
            self.addChildTarget(libsample.SampleAnalysis(sample, samdir,
                                     outfile, sample_sampling_diversity,
                                     size, self.indices, None))
        self.setFollowOnTarget(DiversitySummary(d_dir, self.outdir,
                               self.indices, self.group2samples, self.matched,
                               self.plotfmt, self.pval))

