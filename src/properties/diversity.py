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
import cpickle as pickle

from sonLib.bioio import system

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon


class SampleSamplingDiversity:
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

    def __getitem__(self, name):
        if name not in self.__dict__:
            raise KeyError(("SampleSamplingDiverisity does not have attribute"
                            + "%s\n") % name)
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

def diversity_text(name2sampling, outfile, indices):
    # Rows = Samples, Columns = Diversity_indices
    f = open(outfile, 'w')
    f.write("#Sample\t%s\n" % "\t".join(indices))
    names = sorted(name2sampling.keys())
    for name in names:
        f.write("%s" % name)
        sampling = name2sampling[name]
        for index in indices:
            if index not in sampling:
                f.write("\tNA")
            else:
                f.write("\t%.3f" % sampling[index])
                stdindex = "%s_std" % index
                if stdindex in sampling:
                    std = sampling[stdindex]
                    f.write(" +/- %.3f" % std)
        f.write("\n")
    f.close()

def diversity_latex_tab(f, name2sampling, indices):
    names = sorted(name2sampling.keys())
    for name in names:
        f.write("%s" % name)
        sampling = name2sampling[name]
        for index in indices
            if index not in sampling:
                f.write(" & NA")
            else:
                f.write(" & %.3f" % sampling[index])
                stdindex = "%s_std" % index
                if stdindex in sampling:
                    std = sampling[stdindex]
                    f.write(" \pm %.3f" % std)
        f.write("\\\\\n")
        f.write("\\hline\n")

def diversity_latex(name2sampling, outfile, indices):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Sample"] + indices
    libcommon.tab_header(f, colnames)
    diversity_latex_tab(f, name2sampling, indices)
    caption = "Diversity indices"
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

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
                if index not in sampling:
                    f.write("\tNA")
                else:
                    f.write("\t%.3f" % sampling[index])
        f.write("\n")
    f.close()

def sample_sampling_diversity(sample, size, indices):
    # Sampling sample to "size, and compute the "index" of the subsample
    subsample = libsample.sampling(sample, size)
    sampling = SampleSamplingDiversity()
    counts = [clone.count for clone in subsample.clones]

    import rpy2.rinterface as rinterface
    rinterface.set_initoptions(('rpy2', '--no-save'))
    rinterface.initr()

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
    avrsampling = SampleSamplingDiversity()
    stds = ["%s_std" % i for i in indices]
    for i, index in enumerate(indices):
        vals = [sampling[index] for sampling in samplings]
        avrsampling[index] = np.mean(vals)
        avrsampling[stds[i]] = np.std(vals)
    return avrsampling

def sample_rf_sizes(sample, bin, rf_sizes):
    # Return the list of sampling size for the sample
    sizes = []
    if rf_sizes:
        for s in rf_sizes:
            if s <= sample.size:
                sizes.append(s)
    else:
        if bin > sample.size:
            sizes = [sample.size]
        else:
            sizes = xrange(bin, sample.size, bin)
    return sizes

def check_rarefaction_options(parser, options, args):
    if options.rf_sizes:
        options.rf_sizes = [long(s) for s in options.rf_sizes.split(",")]
    indices = []
    valid_indices = ['numclone', 'simpson', 'invsimpson', 'shannon',
                     'fisher_alpha']
    if options.diversity:
        items = parser.diversity.split(',')
        for item in items:
            if item not in valid_indices:
                raise ValueError("Unknown diversity index: %s" % item)
            indices.append(item)
    options.diversity = indices

def add_rarefaction_options(parser):
    group = OptionGroup(parser, "Rarefaction analyses options")
    group.add_option('--rf_num_sampling', type='int', default=1,
                     help=('Number of samplings performed for each sampling' +
                           'size. Default=%default')
    group.add_option('--bin', dest='bin', type='int', default=50000,
                     help=('Increment of sampling size for rarefaction ' +
                           'analyses. Default=%default, i.e sampling sizes ' +
                           'are 50k, 100k,150k, etc.'))
    group.add_option('--rf_sizes', dest='rf_sizes', 
                     help=('Optional. Comma separated string of sampling ' +
                        'sizes. If specified, will ignore options "bin".')
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
            for file in sizedirpath:
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
            sizedir = os.path.join(self.outdir, size)
            system("mkdir -p %s" % sizedir)
            for i in xrange(self.numsampling):
                outfile = os.path.join(sizedir, "%d.pickle" % i)
                self.addChildTarget(SampleSamplingDiversity(outfile, size,
                                                 self.sample, self.indices)
        self.setFollowOnTarget(AvrSamplingDiversity(self.outdir, self.indices,
                                                    self.avrdir))

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
            sizes = sample_rf_sizes(sample, bin, rf_sizes)
            self.addChildTarget(SampleDiversityRarefaction(sampledir, sample,
                                    sizes, self.options.rf_num_sampling,
                                    self.options.diversity, avrsampledir))
        #self.setFollowOnTarget()

