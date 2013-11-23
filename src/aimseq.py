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


#======== MAIN PIPELINE ========
class Setup(Target):
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
    def __init__(self, name2sam_file, options):
        Target.__init__(self)
        self.name2sam_file = name2sam_file
        self.options = options

    def run(self):
        name2sample = pickle.load(gzip.open(self.name2sam_file, "rb"))
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
        self.setFollowOnTarget()

def add_options(parser):
    group = OptionGroup(parser, "Input arguments")
    group.add_option('-i', '--indir', dest='indir',
                      help=('Input directory containing tab-separated clone '
                            + 'files. Valid formats: MiTCR, Adaptive '
                            + 'Biotechnologies and Sequenta. File names are '
                            + 'used as sample names. (Required argument).'))
    group.add_option('-f', '--format', dest='format', default='mitcr',
                      help= 'Input format. Please choose one of the following:'
                            + ' [mitcr,adaptive,sequenta]. Default=%default'))
    group.add_option('--ext', dest='ext', 
                     help='Input file extension. (Optional)')
    group.add_option('-o', '--outdir', dest='outdir', default='.',
                      help='Output directory. Default=%default')
    group.add_option('-m', '--metadata', dest='metainfo',
                      help=('File containing grouping information. Format:\n'
                            + 'First line:\nmatched=[true/false].Followed by:'
                            + '\n<Group>\\t<sample1,sample2,etc>. One line '
                            + 'group. Note: if matched=true, each group must '
                            + 'have the same number of samples and the order '
                            + 'of the samples implied their matching.'))
    group.add_option('-a', '--analyses', dest='analyses', default='prelim',
                     help=('Types of analyses to perform. Default=%default. '
                           + 'Valid options (comma-separated) are: prelim,'
                           + 'diversity,similarity,clonesize,lendist,geneusage'
                           + ',aausage,trackclone.'))
    group.add_option('--normalize', dest='normalize', action='store_true',
                     default=False, help=('If specified, perform Bioconduct\'s'
                                          + ' metagenomeSeq CSS normalization.'
                                          ))
    group.add_option('--sampling', dest='sampling', type='int',
                     help='Sampling size. Default is using all reads.')
    parser.add_option_group(group)
    incommon.add_filter_options(parser)

def check_options(parser, args, options):
    if options.indir is None:
        raise libcommon.InputError("Please specify input directory.")
    libcommon.check_options_dir(options.indir)
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" % options.outdir)
    libcommon.check_options_dir(options.outdir)
    if options.metainfo:
        libcommon.check_options_file(options.metainfo)

    my_analyses = ['diversity', 'similarity', 'clonesize', 'lendist', 
                   'geneusage', 'aausage', 'trackclone', 'prelim']
    analyses = options.analyses.split(',')
    for a in analyses:
        if a not in my_analyses:
            raise libcommon.InputError("Unknown analysis %s." % a)
    self.analyses = analyses

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

