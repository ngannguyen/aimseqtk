#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Perform rarefaction analyses and compute diversity indices
Plot: xaxis: # sampling size; yaxis: diversity index. One curve/sample
'''

import os
import sys
import numpy as np
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon
import aimseqtk.lib.statcommon as statcommon
import aimseqtk.src.properties.rarefaction_plot as rfplot


class SampleSamplingDiversityStats:
    '''Different indices commpute from samplings of a specific size of
       a sample
    '''
    def __init__(self):
        self.numclone = None
        #R vegan: diversity
        self.simpson = None 
        self.invsimpson = None
        self.shannon = None
        self.fisher_alpha = None

        #Standard deviation
        self.numclone_std = None
        self.simpson_std = None 
        self.invsimpson_std = None
        self.shannon_std = None
        self.fisher_alpha_std = None

    def getitems(self):
        return self.__dict__.keys()
    
    def __getitem__(self, name):
        if name not in self.__dict__:
            return None
            #raise KeyError("SampleSamplingDiverisity does not have attribute %s" % name)
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

def diversity_text(stats, outfile, indices):
    # Rows = Samples, Columns = Diversity_indices
    f = open(outfile, 'w')
    f.write("#Sample\t%s\n" % "\t".join(indices))
    for (name, stat) in stats:
        f.write("%s" % name)
        for index in indices:
            if index not in stat.getitems():
                f.write("\tNA")
            else:
                f.write("\t%.3f" % stat[index])
                stdindex = "%s_std" % index
                if stdindex in stat.getitems():
                    std = stat[stdindex]
                    f.write(" +/- %.3f" % std)
        f.write("\n")
    f.close()

def diversity_latex_tab(f, stats, indices):
    for (name, stat) in stats:
        f.write("%s" % name.replace('_', '\_'))
        for index in indices:
            if index not in stat.getitems():
                f.write(" & NA")
            else:
                f.write(" & %.3f" % stat[index])
                stdindex = "%s_std" % index
                if stdindex in stat.getitems():
                    std = stat[stdindex]
                    f.write(" $\pm$ %.3f" % std)
        f.write("\\\\\n")
        f.write("\\hline\n")

def diversity_latex(stats, outfile, indices):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Sample"] + indices
    colnames = [c.replace('_', '\_') for c in colnames]
    libcommon.tab_header(f, colnames)
    diversity_latex_tab(f, stats, indices)
    caption = "Diversity indices"
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def diversity_table(name2stat, outfile, indices, group2avr={},
                    group2names={}, tex=False, keyindex='numclone'):
    if not name2stat or not indices:
        return
    # Sort the samples by group and/or size
    names = []
    if keyindex not in indices:
        keyindex = indices[0]
    keyfunc = lambda item: item[keyindex]
    if group2avr:
        sortedstats = libcommon.sort_objs_by_group(name2stat, group2names,
                                           True, group2avr, keyfunc=keyfunc)
    else:
        sortedstats = libcommon.sort_dict_by_value(name2stat, keyfunc=keyfunc)
    if tex:
        diversity_latex(sortedstats, outfile, indices)
    else:
        diversity_text(sortedstats, outfile, indices)

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

def sample_sampling_diversity(sample, size, indices):
    # Sampling sample to "size, and compute the "index" of the subsample
    if size is None:  # no sampling, compute diversity of original sample
        subsample = sample
    else:
        subsample = libsample.sampling(sample, size)
    sampling = SampleSamplingDiversityStats()
    counts = [clone.count for clone in subsample.clones]

    #import rpy2.rinterface as rinterface
    #rinterface.set_initoptions(('rpy2', '--no-save', '--no-restore'))
    #rinterface.initr()

    import rpy2.robjects as robjs
    from rpy2.robjects.packages import importr
    vegan = importr("vegan")
    rcounts = robjs.IntVector(counts)
    for index in indices:
        if index == 'numclone':
            sampling[index] = subsample.numclone
        elif index == 'fisher_alpha':
            rfisher = vegan.fisher_alpha(rcounts)
            sampling[index] = rfisher[0]
        else:
            rval = vegan.diversity(rcounts, index)
            sampling[index] = rval[0]
    
    return sampling

def sample_avr_sampling_diversity(samplings, indices):
    avrsampling = SampleSamplingDiversityStats()
    stds = ["%s_std" % i for i in indices]
    for i, index in enumerate(indices):
        vals = [sampling[index] for sampling in samplings]
        avrsampling[index] = np.mean(vals)
        avrsampling[stds[i]] = np.std(vals)
    return avrsampling

def sample_rf_sizes(sample, bin=None, rf_sizes=None):
    # Return the list of sampling size for the sample
    sizes = []
    if rf_sizes:
        for s in rf_sizes:
            if s <= sample.size:
                sizes.append(s)
    elif bin and bin <= sample.size:
        sizes = range(bin, sample.size + 1, bin)
    
    if len(sizes) == 0:
        sizes = [sample.size]
    return sizes

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
    group.add_option('--diversity', dest='diversity', default='numclone',
                     help=('Comma separated list of diversity indices. ' +
                           'Default=%default. Valid options are: [numclone' +
                           ',simpson,invsimpson,shannon,fisher_alpha].'))
    parser.add_option_group(group)

#============ PARALLELIZE =============
class SampleSamplingDiversity(Target):
    def __init__(self, outfile, sample, size, indices):
        Target.__init__(self)
        self.outfile = outfile
        self.sample = sample
        self.size = size
        self.indices = indices

    def run(self):
        sampling = sample_sampling_diversity(self.sample, self.size,
                                             self.indices)
        pickle.dump(sampling, gzip.open(self.outfile, "wb"))

class AvrSamplingDiversity(Target):
    def __init__(self, indir, indices, outdir):
        Target.__init__(self)
        self.indir = indir
        self.indices = indices
        self.outdir = outdir

    def run(self):
        for sizedir in os.listdir(self.indir):
            sizedirpath = os.path.join(self.indir, sizedir)
            samplings = []
            for file in os.listdir(sizedirpath):
                filepath = os.path.join(sizedirpath, file)
                sampling = pickle.load(gzip.open(filepath, "rb"))
                samplings.append(sampling)
            avrsampling = sample_avr_sampling_diversity(samplings, 
                                                        self.indices)
            outfile = os.path.join(self.outdir, "%s.pickle" % sizedir)
            pickle.dump(avrsampling, gzip.open(outfile, "wb"))

class SampleDiversityRarefaction(Target):
    '''Perform rarefaction analyses on the specific sample
    "numsampling" of times
    '''
    def __init__(self, outdir, sample, sizes, numsampling, indices, avrdir):
        Target.__init__(self)
        self.outdir = outdir
        self.sample = sample
        self.sizes = sizes
        self.numsampling = numsampling
        self.indices = indices
        self.avrdir = avrdir

    def run(self):
        for size in self.sizes:
            sizedir = os.path.join(self.outdir, "%d" % size)
            system("mkdir -p %s" % sizedir)
            for i in xrange(self.numsampling):
                outfile = os.path.join(sizedir, "%d.pickle" % i)
                self.addChildTarget(SampleSamplingDiversity(outfile,
                                        self.sample, size, self.indices))
        self.setFollowOnTarget(AvrSamplingDiversity(self.outdir, self.indices,
                                                    self.avrdir))

class DiversityRarefactionSummary(Target):
    # Print summary tables and draw plots
    def __init__(self, indir, name2sample, options):
        Target.__init__(self)
        self.indir = indir   # avrdir: avrdir/sample/size.pickle
        self.name2sample = name2sample
        self.options = options

    def run(self):
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
                rfplot.draw_rarefaction(name2size2sampling, self.name2sample,
                          groups, index, rfplotfile, opts.plotformat, opts.dpi)
        # Summary table for each sampling size:
        for size, name2sampling in size2name2sampling.iteritems():
            g2avr = {}
            if g2s:
                g2avr = libcommon.get_group_avr(name2sampling, g2s)
            txtfile = os.path.join(txtdir, "diversity_%d.txt" % size)
            diversity_table(name2sampling, txtfile, indices, g2avr, g2s)
            texfile = os.path.join(texdir, "diversity_%d.tex" % size)
            diversity_table(name2sampling, texfile, indices, g2avr, g2s, True)

class DiversityRarefaction(Target):
    # diversity_avr/sample/size.pickle
    def __init__(self, name2sample, options):
        Target.__init__(self)
        self.name2sample = name2sample
        self.options = options

    def run(self):
        global_dir = self.getGlobalTempDir()
        outdir = os.path.join(global_dir, "diversity")
        system("mkdir -p %s" % outdir)
        avrdir = os.path.join(global_dir, "diversity_avr")
        system("mkdir -p %s" % avrdir)

        bin = self.options.bin
        rf_sizes = self.options.rf_sizes
        for name, sample in self.name2sample.iteritems():
            sampledir = os.path.join(outdir, name)
            avrsampledir = os.path.join(avrdir, name)
            system("mkdir -p %s" % sampledir)
            system("mkdir -p %s" % avrsampledir)
            sizes = sample_rf_sizes(sample, bin, rf_sizes)
            self.addChildTarget(SampleDiversityRarefaction(sampledir, sample,
                                    sizes, self.options.rf_num_sampling,
                                    self.options.diversity, avrsampledir))
        self.setFollowOnTarget(DiversityRarefactionSummary(avrdir,
                                   self.name2sample, self.options))

class DiversityTtest(Target):
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
        g2stat, g2mean = statcommon.ttest_allpairs(self.g2n, self.name2stat,
                                                   self.matched, self.index)
        picklefile = os.path.join(self.outdir, "%s.pickle" % self.index)
        pickle.dump((g2stat, g2mean), gzip.open(picklefile, 'wb'))
        if self.plotfmt:
            plotfile = os.path.join(plotdir, self.index)
            

class DiversitySummary(Target):
    #  table of diversity indices
    def __init__(self, sampledir, outdir, indices, group2samples, matched,
                                                                 plotfmt=None):
        Target.__init__(self)
        self.sampledir = sampledir
        self.outdir = outdir
        self.indices = indices
        self.group2samples = group2samples
        self.matched = matched
        self.plotfmt = plotfmt

    def run(self):
        name2sampling = {}
        for file in os.listdir(self.sampledir):
            name, ext = os.path.splitext(file)
            filepath = os.path.join(self.sampledir, file)
            sampling = pickle.load(gzip.open(filepath, "rb"))
            name2sampling[name] = sampling
        # Print out summary table
        g2s = self.group2samples
        g2avr = {}
        if g2s:
            g2avr = libcommon.get_group_avr(name2sampling, g2s)
        txtfile = os.path.join(self.outdir, "diversity.txt")
        diversity_table(name2sampling, txtfile, self.indices, g2avr, g2s)
        texfile = os.path.join(self.outdir, "diversity.tex")
        diversity_table(name2sampling, texfile, self.indices, g2avr, g2s, True)

        # For each diversity index, each pair of groups, perform ttest
        # and draw plot
        if g2s and len(g2s) >= 2:
            outdir = os.path.join(self.outdir, "group_comparisons")
            system("mkdir -p %s" % outdir)
            global_dir = self.getGlobalTempDir()
            for index in self.indices:
                self.addChildTarget(DiversityTtest(global_dir, index, g2s,
                            name2sampling, self.matched, self.plotfmt, outdir))
            self.setFollowOnTarget(DiversityTtestSummary(global_dir, outdir,
                                                                 self.indices))

class Diversity(Target):
    '''Calculate diversity indices for input samples. No sampling.
       Table: Rows=Samples; Cols=Indices
       If group2samples is provided: perform statistic test comparing
       pair of groups. If plot is True, draw plot: boxplots of: 
       xaxis: groups, yaxis: diversity indices
    '''
    def __init__(self, samples, outdir, indices, group2samples, matched,
                                                                 plotfmt=None):
        Target.__init__(self)
        self.samples = samples
        self.outdir = outdir
        self.indices = indices
        self.group2samples = group2samples
        self.matched = matched
        self.plotfmt = plotfmt

    def run(self):
        size = None  # no sampling
        global_dir = self.getGlobalTempDir()
        for sample in self.samples:
            outfile = os.path.join(global_dir, "%s.pickle" % sample)
            self.addChildTarget(SampleSamplingDiversity(outfile, sample,
                                                      size, self.indices))
        self.setFollowOnTarget(DiversitySummary(global_dir, self.outdir,
                self.indices, self.group2samples, self.matched, self.plotfmt))

