#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Adaptive IMmune Sequencing ToolKit main pipeline
Inputs:
    1/ Directory containing repertoire TCR/BCR clonotype files
    2/ Sample grouping information
    3/ Analyses wished to performed
Outputs:
    
'''

import os
import sys
import re
import time
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system
from sonLib.bioio import logger

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.drawcommon as drawcommon
import aimseqtk.lib.sample as libsample
import aimseqtk.src.input.inputcommon as incommon
import aimseqtk.src.properties.repsize as repsize
import aimseqtk.src.properties.diversity as diversity
import aimseqtk.src.normalize as normalize


#======== MAIN PIPELINE ========
class Setup(Target):
    '''Sets up AIMSEQTK pipeline, starting with parsing input clone files
       If "pickle" is the input format, will assume that "prelim" steps
       were done, and proceed to "Preprocess"
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        if self.options.format == 'pickle':
            self.addChildTarget(Preprocess(sampledir, self.options))
        else:
            # Read input files:
            global_dir = self.getGlobalTempDir()
            sampledir = os.path.join(global_dir, "samples")
            system("mkdir -p %s" % sampledir)
            self.addChildTarget(incommon.ReadCloneFiles(sampledir,
                                self.options.indir, self.options.format,
                                self.options.ext))
            self.setFollowOnTarget(Filter(sampledir, self.options))

class Filter(Target):
    '''Filter each sample by size (min and/or max count and/or freq)
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        opts = self.options
        name2group = {}
        name2color = {}
        name2marker = {}
        if opts.group2samples:
            name2group = libcommon.get_val2key_1to1(opts.group2samples)
        if opts.makeplots:
            names = [os.path.splitext(file)[0] for file in 
                                                os.listdir(self.sampledir)]
            assert names is not None
            name2color = drawcommon.get_name2color_wtgroup(names, name2group,
                                                   opts.group2samples)
            name2marker = drawcommon.get_name2marker_wtgroup(names,
                                                   opts.group2samples)
        global_dir = self.getGlobalTempDir()
        filter_dir = os.path.join(global_dir, "filterBySize")
        system("mkdir -p %s" % filter_dir)
        for file in os.listdir(self.sampledir):
            samplefile = os.path.join(self.sampledir, file)
            sample = pickle.load(gzip.open(samplefile, "rb"))
            libsample.set_sample_group(sample, name2group)
            libsample.set_sample_color(sample, name2color)
            libsample.set_sample_marker(sample, name2marker)
            outfile = os.path.join(filter_dir, "%s.pickle" % sample.name)
            self.addChildTarget(libsample.FilterBySize(outfile, sample,
                   opts.mincount, opts.maxcount, 
                   opts.minfreq, opts.maxfreq, True))
        self.setFollowOnTarget(GetProductive(filter_dir, opts))

class GetProductive(Target):
    '''Split each sample into productive and non-productive clones
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        global_dir = self.getGlobalTempDir()
        outdir = os.path.join(global_dir, "filterByStatus")
        pdir = os.path.join(outdir, "productive")
        system("mkdir -p %s" % pdir)
        npdir = os.path.join(outdir, "non_productive")
        system("mkdir -p %s" % npdir)

        for file in os.listdir(self.indir):
            filepath = os.path.join(self.indir, file)
            sample = pickle.load(gzip.open(filepath, "rb"))
            self.addChildTarget(libsample.FilterByStatus(pdir, npdir, sample,
                                                         True, self.options))
        self.setFollowOnTarget(Preprocess(pdir, self.options))

class Preprocess(Target):
    '''If at the "prelim" stage, perform rarefaction analyses
       If "sampling" is requested, down-sampling samples, then 
       if "normalize" is requested, normalize samples, then
       go to the beginning of analyses
    '''
    def __init__(self, sampledir, options):
        self.sampledir = sampledir
        self.options = options

    def run(self):
        opts = self.options
        if opts.analyses == ['prelim']:
            self.addChildTarget(RepSize(self.sampledir, opts))
        elif opts.sampling:
            self.addChildTarget(DownSampling(self.sampledir, opts))
        elif opts.normalize:
            self.addChildTarget(Normalize(self.sampledir, opts))
        else:
            self.addChildTarget(Analyses(self.sampledir, opts))
        #self.setFollowOnTarget()

class RepSize(Target):
    '''Summarize samples' #clones & # reads 
       Rarefaction analyses
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = indir
        self.options = options

    def run(self):
        name2sample = {}
        for file in os.listdir(self.sampledir):
            filepath = os.path.join(self.sampledir, file)
            prd_sam = pickle.load(gzip.open(filepath, 'rb'))
            name2sample[prd_sam.name] = prd_sam

        # Get summary of samples' sizes:
        group2samples = self.options.group2samples
        #group2avr = repsize.get_group_avr(name2sample, group2samples)
        group2avr = libcommon.get_group_avr(name2sample, group2samples)
        txtfile = os.path.join(self.options.outdir, "clonesize.txt")
        repsize.repsize_table(name2sample, txtfile, group2avr, group2samples)
        texfile = os.path.join(self.options.outdir, "clonesize.tex")
        repsize.repsize_table(name2sample, texfile, group2avr, group2samples,
                              True)
        self.addChildTarget(diversity.DiversityRarefaction(name2sample,
                                                           self.options))

class DownSampling(Target):
    '''Down-sampling each sample to a specific size, then continue to
       normalization step if requestedd, otherwise move to analyses
    '''
    def __init__(self, sampledir, options):
        Target.__init(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        opts = self.options
        global_dir = self.getGlobalTempDir()
        for file in os.listdir(self.sampledir):
            filepath = os.path.join(self.sampledir, file)
            sample = pickle.load(gzip.open(filepath, 'rb'))
            self.addChildTarget(libsample.Sampling(sample, opts.sampling,
                                                   global_dir))
        if opts.normalize:
            self.setFollowOnTarget(Normalize(global_dir, opts))
        else:
            self.setFollowOnTarget(Analyses(global_dir, opts))

class Normalize(Target):
    '''Normalize clones' sizes and continue to analyses
    '''
    def __init__(self, sampledir, options):
       Target.__init__(self)
       self.sampledir = sampledir
       self.options = options

    def run(self):
        opts = self.options
        global_dir = self.getGlobalTempDir()
        samples = libcommon.load_pickledir(self.sampledir)
        sizetype = 'count'
        norm_matrix = normalize.normalize_MRexp(samples, sizetype)
        samples = normalize.matrix_to_normcount(norm_matrix, samples)
        for s in samples:
            outfile = os.path.join(global_dir, "%s.pickle" % s.name)
            pickle.dump(s, gzip.open(outfile, "wb"))
        self.setFollowOnTarget(Analyses(global_dir, opts))

class Analyses(Target):
    '''
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        opts = self.options
        samples = libcommon.load_pickledir(self.sampledir)
        plotfmt = None
        if opts.makeplots:
            plotfmt = opts.plotformat
        if 'diversity' in opts.analyses:
            diversitydir = analysis_outdir("diversity", opts.outdir)
            self.addChildTarget(diversity.Diversity(samples, diversitydir,
                                       opts.diversity, opts.group2samples,
                                       opts.matched, plotfmt))
        if 'similarity' in opts.analyses:
            similaritydir = analysis_outdir("similarity", opts.outdir)
            self.addChildTarget(similarity.Similarity(samples, similaritydir,
                                                      opts))
        if 'clonesize' in opts.analyses:
            clonesizedir = analysis_outdir("clonesize", opts.outdir)
            self.addChildTarget(clonesize.CloneSize(samples, clonesizedir,
                                                                      opts))
        if 'lendist' in opts.analyses:
            lendistdir = analysis_outdir("lendist", opts.outdir)
            self.addChildTarget(lendist.LenDist(samples, lendistdir, opts))
        if 'geneusage' in opts.analyses:
            geneusagedir = analysis_outdir("geneusage", opts.outdir)
            self.addChildTarget(geneusage.GeneUsage(samples, geneusagedir,
                                                                     opts))
        if 'aausage' in opts.analyses:
            aausagedir = analysis_outdir("aausage", opts.outdir)
            self.addChildTarget(aausage.AaUsage(samples, aausagedir, opts))
        if 'trackclone' in opts.analyses:
            trackclonedir = analysis_outdir("trackclone", opts.outdir)
            self.addChildTarget(trackclone.TrackClone(samples, trackclonedir,
                                                                       opts))

def analysis_outdir(analysis, outdir):
    analysis_dir = os.path.join(outdir, analysis)
    system("mkdir -p %s" % analysis_dir)
    return analysis_dir

def add_options(parser):
    incommon.add_input_options(parser)
    incommon.add_filter_options(parser)
    diversity.add_rarefaction_options(parser)

def check_options(parser, options):
    incommon.check_input_options(parser, options)
    diversity.check_rarefaction_options(parser, options)

def main():
    parser = libcommon.init_options()
    add_options(parser)
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    check_options(parser, options)

    i = Stack(Setup(options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__":
    from aimseqtk.src.aimseq import *
    main()

