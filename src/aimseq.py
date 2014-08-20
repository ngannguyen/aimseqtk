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
import aimseqtk.src.properties.similarity as similarity
import aimseqtk.src.properties.clonesize as clonesize
import aimseqtk.src.normalize.normalize as normalize
import aimseqtk.src.cdr3len.cdr3len as lendist
import aimseqtk.src.geneusage.geneusage as geneusage
import aimseqtk.src.overlap.overlap as overlap
import aimseqtk.src.overlap.trackclone as trackclone
import aimseqtk.src.recomb.recomb_model as recomb_model


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
        self.logToMaster("Setting up...")
        opts = self.options
        indir = opts.indir
        if opts.regroup:  # regrouping the samples
            indir = opts.regroup
            samples_regroup(opts.indir, indir, opts.group2samples, opts.groups, opts.matched)

        if opts.format == 'pickle' or \
           (opts.format == 'db' and 'model' in opts.analyses):
            self.addChildTarget(MakeDbSamples(indir, opts))
        elif opts.format == 'db':
            self.addChildTarget(Preprocess(indir, opts))
        else:
            # Read input files:
            global_dir = self.getGlobalTempDir()
            sampledir = os.path.join(global_dir, "samples")
            system("mkdir -p %s" % sampledir)
            self.addChildTarget(incommon.ReadCloneFiles(sampledir,
                                indir, opts.format, opts.ext))
            self.setFollowOnTarget(Filter(sampledir, opts))

class Filter(Target):
    '''Filter each sample by size (min_ and/or max_ count and/or freq)
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        opts = self.options
        global_dir = self.getGlobalTempDir()
        filter_dir = os.path.join(global_dir, "filtered")
        system("mkdir -p %s" % filter_dir)
        for sam in os.listdir(self.sampledir):  # each sample
            samdir = os.path.join(self.sampledir, sam)  
            if not os.path.isdir(samdir):
                raise incommon.ReadSampleError("%s is not a dir" % samdir)
            for file in os.listdir(samdir):  # each batch
                samfile = os.path.join(samdir, file)
                self.addChildTarget(libsample.FilterSample(filter_dir, sam,
                                                           samfile, opts))
        self.setFollowOnTarget(MakeDbSamples(filter_dir, opts))

class MakeDbSamples(Target):
    '''Set up jobs to load each sample into one sqlite3 db file,
    which contains:
    1/ a sampleinfo table,
    2/ clone tables, each represent a V-J combination
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        self.logToMaster("MakeDbSamples\n")
        opts = self.options

        make_db = False
        if 'db' in opts.analyses:
            make_db = True
        #make_db = True
        get_model = False
        if 'model' in opts.analyses:
            get_model = True

        dbdir = os.path.join(opts.outdir, "sample_db")
        if make_db:
            system("mkdir -p %s" % dbdir)
        modeldir = os.path.join(opts.outdir, "recomb_model")
        if get_model:
            system("mkdir -p %s" % modeldir)

        #status = ['productive']
        status = ['productive', 'non_productive']
        for s in status:
            indir = os.path.join(self.sampledir, s)
            names = os.listdir(indir)
            outdir = os.path.join(dbdir, s)
            n2g, n2c, n2m = samples_group_info(names, opts.group2samples, opts.groups, opts.matched)
            if make_db:
                system("mkdir -p %s" % outdir)
            
            for sam in names:
                samdir = os.path.join(indir, sam)
                assert os.path.isdir(samdir)
                g, c, m = sample_group_info(sam, n2g, n2c, n2m)
                if make_db:
                    sam_dbdir = os.path.join(outdir, sam)
                    system("mkdir -p %s" % sam_dbdir)
                    self.addChildTarget(libsample.MakeDbSample(samdir, sam_dbdir, opts,
                                                           g, c, m))
                if get_model:
                    sam_modeldir = os.path.join(modeldir, s, sam)
                    system("mkdir -p %s" % sam_modeldir)
                    self.addChildTarget(recomb_model.SampleRecombModel(samdir,
                                                                 sam_modeldir))
            if opts.samout:  # write samples if requested
                samoutdir = os.path.join(opts.outdir, "samples", s)
                system("mkdir -p %s" % samoutdir)
                self.addChildTarget(libsample.WriteSamples(indir,
                                                    samoutdir, opts.samout))
        if (opts.analyses != ['db'] and
            set(opts.analyses) != set(['db', 'model']) and
            opts.analyses != ['model']):
            pdir = os.path.join(dbdir, "productive")
            self.setFollowOnTarget(Preprocess(pdir, opts))

class Preprocess(Target):
    '''If at the "prelim" stage, perform rarefaction analyses
       If "sampling" is requested, down-sampling samples, then 
       if "normalize" is requested, normalize samples, then
       go to the beginning of analyses
       sampledir/
          sample1/
             sample1
             v1j1
             viji
             ...
             vnjn
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        self.logToMaster("Preprocess\n")
        opts = self.options
        #if opts.analyses == ['prelim']:
        if 'prelim' in opts.analyses:
            self.addChildTarget(RepSize(self.sampledir, opts))
        elif opts.sampling or opts.sampling_uniq:
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
        self.sampledir = sampledir
        self.options = options

    def run(self):
        self.logToMaster("RepSize\n")
        stime = time.time()
        name2sample = {}
        for sam in os.listdir(self.sampledir):
            filepath = os.path.join(self.sampledir, sam, sam)
            sample = pickle.load(gzip.open(filepath, 'rb'))
            name2sample[sam] = sample
        logger.info("RepSize, done loading %d samples in %.4f s." %
                    (len(name2sample), (time.time() - stime)))
        stime = time.time()

        # Get summary of samples' sizes:
        group2samples = self.options.group2samples
        group2avr = libcommon.get_group_avr(name2sample, group2samples)
        logger.info("RepSize, done computing group_avr in %.4f s." %
                    (time.time() - stime))
        
        txtfile = os.path.join(self.options.outdir, "clonesize.txt")
        repsize.repsize_table(name2sample, txtfile, group2avr, group2samples)
        texfile = os.path.join(self.options.outdir, "clonesize.tex")
        repsize.repsize_table(name2sample, texfile, group2avr, group2samples,
                              True)
        self.addChildTarget(diversity.DiversityRarefaction(self.sampledir,
                                                           self.options))

class DownSampling(Target):
    '''Down-sampling each sample to a specific size, then continue to
       normalization step if requestedd, otherwise move to analyses
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        self.logToMaster("DownSampling\n")
        opts = self.options
        global_dir = self.getGlobalTempDir()
        #sampling_dir = os.path.join(global_dir, "down_sampling")
        sampling_dir = os.path.join(opts.outdir, "down_sampling")
        system("mkdir -p %s" % sampling_dir)

        for sam in os.listdir(self.sampledir):
            samdir = os.path.join(self.sampledir, sam)
            sample = pickle.load(gzip.open(os.path.join(samdir, sam), "rb"))
            out_samdir = os.path.join(sampling_dir, sam) 
            system("mkdir -p %s" % out_samdir)
            if opts.sampling_uniq:  # sampling uniq clones
                self.addChildTarget(libsample.SampleAnalysis0(sample, samdir,
                                              out_samdir, libsample.sampling,
                                              opts.sampling_uniq, 'uniq'))
            elif opts.sampling_top:  # sampling reads, then report top clones
                self.addChildTarget(libsample.SampleAnalysis0(sample, samdir,
                                out_samdir, libsample.sampling, opts.sampling,
                                "top", opts.sampling_top))
            else:  # sampling reads
                self.addChildTarget(libsample.SampleAnalysis0(sample, samdir,
                                out_samdir, libsample.sampling, opts.sampling))
        if opts.normalize:
            self.setFollowOnTarget(Normalize(sampling_dir, opts))
        else:
            self.setFollowOnTarget(Analyses(sampling_dir, opts))

class Normalize(Target):
    '''Normalize clones' sizes and continue to analyses
    '''
    def __init__(self, sampledir, options):
       Target.__init__(self)
       self.sampledir = sampledir
       self.options = options

    def run(self):
        self.logToMaster("Normalize\n")
        opts = self.options
        global_dir = self.getGlobalTempDir()
        #norm_dir = os.path.join(global_dir, "normalized")
        norm_dir = os.path.join(opts.outdir, "normalized")
        system("mkdir -p %s" % norm_dir)

        #samples = libcommon.load_pickledir(self.sampledir)
        sizetype = 'count'
        self.addChildTarget(normalize.NormalizeMRexp(self.sampledir, norm_dir,
                                                     sizetype))
        self.setFollowOnTarget(Analyses(norm_dir, opts))

class Analyses(Target):
    '''
    '''
    def __init__(self, sampledir, options):
        Target.__init__(self)
        self.sampledir = sampledir
        self.options = options

    def run(self):
        self.logToMaster("Analyses\n")
        opts = self.options
        samdir = self.sampledir
        #samples = libcommon.load_pickledir(self.sampledir)
        plotfmt = None
        if opts.makeplots:
            plotfmt = opts.plotformat
        if 'diversity' in opts.analyses:
            diverdir = analysis_outdir("diversity", opts.outdir)
            self.addChildTarget(diversity.Diversity(samdir, diverdir,
                                       opts.diversity, opts.group2samples,
                                       opts.matched, plotfmt, opts.pval))
        if 'similarity' in opts.analyses:
            simidir = analysis_outdir("similarity", opts.outdir)
            self.addChildTarget(similarity.Similarity(samdir, simidir, opts))
        if 'clonesize' in opts.analyses:
            csdir = analysis_outdir("clonesize", opts.outdir)
            self.addChildTarget(clonesize.CloneSize(samdir, csdir, opts))
        if 'lendist' in opts.analyses:
            lddir = analysis_outdir("lendist", opts.outdir)
            self.addChildTarget(lendist.LenDist(samdir, lddir, opts))
        if 'geneusage' in opts.analyses:
            gudir = analysis_outdir("geneusage", opts.outdir)
            self.addChildTarget(geneusage.GeneUsage(samdir, gudir, opts))
        #if 'aausage' in opts.analyses:
        #    audir = analysis_outdir("aausage", opts.outdir)
        #    self.addChildTarget(aausage.AaUsage(samples, audir, opts))
        if 'overlap' in opts.analyses:
            odir = analysis_outdir("overlap", opts.outdir)
            self.addChildTarget(overlap.Overlap(samdir, odir, opts))
        if 'trackclone' in opts.analyses:
            tdir = analysis_outdir("track_top_clones", opts.outdir)
            self.addChildTarget(trackclone.TrackTopClones(samdir, tdir, opts))
        if opts.clones:
            cdir = analysis_outdir("track_clones", opts.outdir)
            self.addChildTarget(trackclone.TrackClones(samdir, cdir, opts,
                                                                  opts.clones))

def sample_group_info(name, n2g, n2c, n2m):
    group = None
    color = (0, 0, 0)
    marker = '.'
    if name in n2g:
        group = n2g[name]
    if name in n2c:
        color = n2c[name]
    if name in n2m:
        marker = n2m[name]
    return group, color, marker

def samples_group_info(names, g2names, groups, matched):
    n2group = {}
    n2color = {}
    n2marker = {}
    if g2names:
        n2group = libcommon.get_val2key_1to1(g2names)
        n2color = drawcommon.get_name2color_wtgroup(names, n2group, g2names, groups, matched)
        n2marker = drawcommon.get_name2marker_wtgroup(names, g2names)
    return n2group, n2color, n2marker

def samples_regroup(indir, outdir, g2s, groups, matched):
    names = os.listdir(indir)
    n2g, n2c, n2m = samples_group_info(names, g2s, groups, matched)
    system('mkdir -p %s' % outdir)
    for sam in os.listdir(indir):
        g, c, m = sample_group_info(sam, n2g, n2c, n2m)
        samdir = os.path.join(indir, sam)
        outsamdir = os.path.join(outdir, sam)
        system('mkdir -p %s' % outsamdir)
        for file in os.listdir(samdir):
            filepath = os.path.abspath(os.path.join(samdir, file))
            outfile = os.path.join(outsamdir, file)
            if file == sam:
                sample = pickle.load(gzip.open(filepath, 'rb'))
                sample.group = g
                sample.color = c
                sample.marker = m
                pickle.dump(sample, gzip.open(outfile, 'wb'))
            else:
                system('ln -s %s %s' % (filepath, outfile))

def analysis_outdir(analysis, outdir):
    analysis_dir = os.path.join(outdir, analysis)
    system("mkdir -p %s" % analysis_dir)
    return analysis_dir

def add_options(parser):
    incommon.add_input_options(parser)
    incommon.add_filter_options(parser)
    diversity.add_rarefaction_options(parser)
    similarity.add_similarity_options(parser)
    clonesize.add_clonesize_options(parser)
    overlap.add_overlap_options(parser)
    trackclone.add_track_clone_options(parser)

def check_options(parser, options):
    incommon.check_input_options(parser, options)
    diversity.check_rarefaction_options(parser, options)
    similarity.check_similarity_options(options)
    #clonesize.check_clonesize_options(parser, options)
    trackclone.check_track_clone_options(parser, options)

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

