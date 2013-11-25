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
import cpickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.sample as libsample
import aimseqtk.src.input.inputcommon as incommon
import aimseqtk.src.properties.repsize as repsize
import aimseqtk.src.properties.diversity as diversity


#======== MAIN PIPELINE ========
class Setup(Target):
    '''Sets up AIMSEQTK pipeline, starting with parsing input clone files
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        #Read input files:
        global_dir = self.getGlobalTempDir()
        nam2sam_file = os.path.join(global_dir, "name2sample.pickle")
        self.addChildTarget(incommon.ReadCloneFiles(nam2sam_file,
                            self.options.indir, self.options.format,
                            self.options.ext))
        self.setFollowOnTarget(Filter(name2sam_file, self.options))

class Filter(Target):
    '''Filter each sample by size (min and/or max count and/or freq)
    '''
    def __init__(self, name2sam_file, options):
        Target.__init__(self)
        self.name2sam_file = name2sam_file
        self.options = options

    def run(self):
        name2sample = pickle.load(gzip.open(self.name2sam_file, "rb"))
        if self.options.group2samples:
            incommon.set_sample_group(name2sample, self.options.group2samples)
        if self.options.makeplots:
            incommon.set_sample_color(name2sample, self.options.group2samples)
            incommon.set_sample_marker(name2sample, self.options.group2samples)
        
        global_dir = self.getGlobalTempDir()
        filter_dir = os.path.join(global_dir, "filterBySize")
        system("mkdir -p %s" % filter_dir)
        for name, sample in name2sample.iteritems():
            outfile = os.path.join(filter_dir, "%s.pickle" % name)
            self.addChildTarget(libsample.FilterBySize(outfile, sample,
                   self.options.mincount, self.options.maxcount, 
                   self.options.minfreq, self.options.maxfreq, True))
        self.setFollowOnTarget(GetProductive(filter_dir, self.options))

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
        system("mkdir -p %s" % outdir)

        files = os.listdir(self.indir)
        for file in files:
            filepath = os.path.join(self.indir, file)
            sample = pickle.load(gzip.open(filepath, "rb"))
            outfile = os.path.join(outdir, file)
            self.addChildTarget(libsample.FilterByStatus(outfile, sample,
                                                         True))
        self.setFollowOnTarget(SamplingAndAnalyses())

class RepSize(Target):
    '''Summarize samples' #clones & # reads 
       Rarefaction analyses
    '''
    def __init__(self, indir, options):
        Target.__init__(self)
        self.indir = indir
        self.options = options

    def run(self):
        name2sample = {}
        for file in os.listdir(self.indir):
            filepath = os.path.join(self.indir, file)
            (prd_sam, nonprd_sam) = pickle.load(gzip.open(file, 'rb'))
            name2sample[prd_sam.name] = prd_sam

        # Get summary of samples' sizes:
        group2samples = self.options.group2samples
        group2avr = repsize.get_group_avr(name2sample, group2samples)
        txtfile = os.path.join(self.options.outdir, "clonesize.txt")
        repsize.repsize_table(name2sample, txtfile, group2avr, group2samples)
        texfile = os.path.join(self.options.outdir, "clonesize.tex")
        repsize.repsize_table(name2sample, txtfile, group2avr, group2samples,
                              True)
        self.addChildTarget(diversity.DiversityRarefaction(name2sample,
                                                           self.options)) 

def add_options(parser):
    incommon.add_input_options(parser)
    incommon.add_filter_options(parser)

def check_options(parser, args, options):
    incommon.check_input_options(parser, args, options)

def main():
    parser = libcommon.init_options()
    add_options()
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    check_options(parser, args, options)

    i = Stack(Setup(options)).startJobTree()
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__":
    from aimseqtk.src.aimseq import *
    main()

