#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys
import time
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system
from sonLib.bioio import logger

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.drawcommon as drawcommon
import aimseqtk.lib.sample as libsample
import aimseqtk.lib.clone as libclone
import aimseqtk.src.input.mitcr as mitcr
import aimseqtk.src.input.adaptive as adaptive
import aimseqtk.src.input.sequenta as sequenta


class FormatError(Exception):
    pass

class ReadSampleError(Exception):
    pass

class UnrecognizedFormatError(Exception):
    pass

def split_file(file, lineperfile, outbase):
    headerfile = "%s_header.txt" % outbase
    system("head -n 1 %s > %s" % (file, headerfile))
    contentfile = "%s_content.txt" % outbase
    system("tail -n +2 %s > %s" % (file, contentfile))
    system("splitFile -head=%s %s %d %s" %
                            (headerfile, contentfile, lineperfile, outbase))
    system("rm -f %s %s" % (headerfile, contentfile))

def read_clone_file(file, parsefunc=mitcr.mitcr_parseline):
    #samplename = os.path.splitext(os.path.basename(file))[0]
    clones = []
    
    f = open(file, 'r')
    headerline = f.readline().lstrip('#')
    index2col = libcommon.get_index2item(headerline)
    for line in f:
        clone = parsefunc(line, index2col)
        if clone is not None:
            clones.append(clone)
    f.close()

    if not clones:
        raise FormatError("File %s has zero clone. Please check the header \
                           line for appropriate column names." % file)

    return clones

def read_clone_files(indir, parsefunc=mitcr.mitcr_parseline, ext=None):
    # "func" is the converting function, different input formats
    name2sample = {}
    for file in os.listdir(indir):
        # Check for the correct file extension if ext is specified
        if ext is not None:
            basename, extension = os.path.splitext(file)
            if not extension.endswith(ext):
                continue

        clones = read_clone_file(os.path.join(indir, file), parsefunc)
        samplename = os.path.splitext(file)[0]
        sample = libsample.Sample(samplename, clones)
        if sample.name in name2sample:
            sys.stderr.write("Warning: Repetitive sample %s\n" % sample.name)
            name2sample[sample.name].addclones(sample.clones)
        else:
            name2sample[sample.name] = sample
    
    return name2sample
        
def process_clone_inputs(indir, format="mitcr", ext=None):
    # format has to be one of the following values: mitcr, adaptive, sequenta
    # return name2sample (key = sample.name, val = Sample())
    format2func = {"mitcr": mitcr.mitcr_parseline,
                   "adaptive": adaptive.adaptive_parseline,
                   "sequenta": sequenta.sequenta_parseline,
                   "aimseqtk": libclone.clone_parseline}
    if format not in format2func:
        types = ",".join(format2func.keys())
        raise UnrecognizedFormatError("Format %s is not recognized. Please \
                    choose one of these: %s\n" % (format, types))
    else:
        func = format2func[format]
        return read_clone_files(indir, func, ext)
    return None

def read_group_info(file):
    # The input file has the following format:
    # matched=[true/false]
    # group1\tsample1_1,sample1_2,sample1_3,etc
    # group2\tsample2_1,sample2_2,etc
    #
    # "matched" is useful in time series or case/control analyses
    # Note: if matched=true, each group must have the same # of samples
    # And the order of the samples implied their matching, 
    # e.g sample1_1 matches sample2_1
    # If matched=false, the # of samples in group can vary and 
    # the ordering is not important
    f = open(file, 'r')
    group2samples = {}
    groups = []
    
    firstline = f.readline()
    if firstline is None:
        raise FormatError("Empty group_info file %s\n" % file)
    firstitems = firstline.strip('\n').split('=')
    if (len(firstitems) != 2 or firstitems[0] != 'matched' or 
               firstitems[1].lower() not in ['true', 'false']):
        raise FormatError("Wrong group_info file format; 1st line must be:\n"
                          + "matched=[true/false]\n. 1st line read was:\n"
                          + "%s\n" % firstline)
    matched = False
    if firstitems[1].lower() == 'true':
        matched = True

    for line in f:
        items = line.strip('\n').split()
        if len(items) != 2:
            raise FormatError("Wrong group_info file format; except for the "
                              + "1st line, each other line must have 2 fields"
                              + ". Line %s has %d field(s)\n" % 
                              (line, len(items)))
        group = items[0]
        samples = items[1].split(',')
        if not samples:
            sys.stderr.write("Group %s has 0 sample.\n" % group)
        if group in group2samples:
            sys.stderr.write("Warning: group_info file %s has groups with the "
                             + "same name: %s -- merge them into 1 group.\n" % 
                             (file, group))
            group2samples[group].extend(samples)
        else:
            group2samples[group] = samples
        groups.append(group)

    if len(group2samples) == 0:
        raise FormatError("Wrong group_info format: No group was specified.\n")
    
    # Check the constraints of matched samples:
    if matched:
        numsams = [len(samples) for samples in group2samples.values()]
        nonredundant = set(numsams)
        if len(nonredundant) != 1:
            raise FormatError("Wrong group_info format. 'matched=true' "
                              + "requires that all group must have the same"
                              + "number of samples.\n")

    f.close()
    return groups, group2samples, matched

#============ OPTIONS ==============
def check_input_options(parser, options):
    if options.indir is None:
        raise libcommon.InputError("Please specify input directory.")
    libcommon.check_options_dir(options.indir)
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" % options.outdir)
    libcommon.check_options_dir(options.outdir)
    
    if not options.jobTree:
        options.jobTree = os.path.join(options.outdir, "jobTree")
    
    options.group2samples = None
    options.matched = None
    if options.metainfo:
        libcommon.check_options_file(options.metainfo)
        groups, group2samples, matched = read_group_info(options.metainfo)
        options.group2samples = group2samples
        options.matched = matched
        options.groups = groups

    my_analyses = ['diversity', 'similarity', 'clonesize', 'lendist', 
                   'geneusage', 'aausage', 'overlap', 'trackclone',
                   'prelim', 'db', 'model']
    analyses = []
    if options.analyses:
        analyses = options.analyses.split(',')
    if 'all' in analyses:
        analyses = my_analyses[:-3]
        #if not options.matched:
        #    analyses.remove('trackclone')
    else:
        for a in analyses:
            if a not in my_analyses:
                raise libcommon.InputError("Unknown analysis %s." % a)
        #if 'trackclone' in analyses and not options.matched:
        #    raise libcommon.InputError(("Matched samples are required if " +
        #                                "'trackclone' analysis is specified." +
        #                                " Please modify the metadata file."))
    options.analyses = analyses
    if options.samout:
        samouts = options.samout.split(',')
        for samout in samouts:
            if samout not in ['fasta', 'txt', 'pickle', 'all']:
                raise libcommon.InputError("Unknown sample out format %s." %
                                            samout)
        if 'all' in samouts:
            options.samout = ['fasta', 'txt', 'pickle']
        else:
            options.samout = samouts
        if 'prelim' in analyses and 'pickle' not in options.samout:
            options.samout.append('pickle')
        samoutdir = os.path.join(options.outdir, "samples")
        for format in options.samout:  # e.g: outdir/samples/productive/pickle
            dir1 = os.path.join(samoutdir, "productive", format)
            system("mkdir -p %s" % dir1)
            dir2 = os.path.join(samoutdir, "non_productive", format)
            system("mkdir -p %s" % dir2)
    drawcommon.check_plot_options(parser, options)

def add_input_options(parser):
    group = OptionGroup(parser, "Input arguments")
    group.add_option('-i', '--indir', dest='indir',
                      help=('Input directory containing tab-separated clone '
                            + 'files. Valid formats: MiTCR, Adaptive '
                            + 'Biotechnologies and Sequenta. File names are '
                            + 'used as sample names. (Required argument).'))
    group.add_option('-f', '--format', dest='format', default='mitcr',
                      help=('Input format. Please choose one of the following:'
                            + ' [mitcr,adaptive,sequenta,aimseqtk,pickle,db]. '
                            + 'Default=%default'))
    group.add_option('--regroup', dest='regroup', default=None,
                      help=('If input format is "db" and need to regroup' +
                            ' the samples. This is the new dbdir.'))
    group.add_option('--ext', dest='ext', 
                     help='Input file extension. (Optional)')
    group.add_option('-o', '--outdir', dest='outdir', default='.',
                      help='Output directory. Default=%default')
    group.add_option('--sample_out_format', dest='samout', help=('Optional: ' +
                     '[fasta,txt,pickle,both]. If specified, will print out ' +
                     'processed (after size and status filtering) samples in' +
                     'the chosen format to outdir/samples. Default: do not ' +
                     'print out samples.'))
    group.add_option('-m', '--metadata', dest='metainfo',
                      help=('File containing grouping information. Format:\n'
                            + 'First line:\nmatched=[true/false].Followed by:'
                            + '\n<Group>\\t<sample1,sample2,etc>. One line '
                            + 'group. Note: if matched=true, each group must '
                            + 'have the same number of samples and the order '
                            + 'of the samples implied their matching.'))
    group.add_option('-a', '--analyses', dest='analyses', default='db',
                     help=('Types of analyses to perform. Default=%default. '
                           + 'Valid options (comma-separated) are: db,prelim,'
                           + 'diversity,similarity,clonesize,lendist,geneusage'
                           + ',aausage,overlap,trackclone,model.'))
    group.add_option('--normalize', dest='normalize', action='store_true',
                     default=False, help=('If specified, perform Bioconduct\'s'
                                          + ' metagenomeSeq CSS normalization.'
                                          ))
    group.add_option('--sampling', dest='sampling', type='long',
                     help='Sampling size. Default is using all reads.')
    group.add_option('--sampling_uniq', dest='sampling_uniq', type='long',
                     help='Sampling uniq clones. Default is using all clones.')
    group.add_option('--sampling_top', dest='sampling_top', type='long',
                     help=('Sampling, then only return the top clones. ' +
                           'Default is using all clones.'))
    group.add_option('--pval', dest='pval', type='float', default=0.05,
                     help='pvalue cutoff')
    parser.add_option_group(group)
    drawcommon.add_plot_options(parser)

def add_filter_options(parser):
    group = OptionGroup(parser, "Repertoire filtering options")
    group.add_option('--mincount', dest='mincount', type='int', default=1,
                      help=('Minimum read count. Clones with smaller counts '
                            + 'will be filtered out. Default=%default'))
    group.add_option('--maxcount', dest='maxcount', type='int', default=-1,
                      help=('Maximum read count. Clones with larger counts '
                            + 'will be filtered out. Default=%default (None).')
                    )
    group.add_option('--minfreq', dest='minfreq', type='float', default=0.0,
                      help=('Minimum read freq. Clones with smaller freqs '
                            + 'will be filtered out. Default=%default'))
    group.add_option('--maxfreq', dest='maxfreq', type='float', default=-1.0,
                      help=('Maximum read freq. Clones with larger freqs '
                            + 'will be filtered out. Default=%default (None).')
                    )
    parser.add_option_group(group)

#============= Parallelizing =========
class ReadCloneFileNoSplit(Target):
    def __init__(self, infile, outfile, parsefunc):
        Target.__init__(self)
        self.infile = infile
        self.outfile = outfile
        self.parsefunc = parsefunc

    def run(self):
        starttime = time.time()
        clones = read_clone_file(self.infile, self.parsefunc)
        pickle.dump(clones, gzip.open(self.outfile, "wb"))
        mytime = time.time() - starttime
        msg = "Done reading sample file %s in %.4f s." % (self.infile, mytime)
        self.logToMaster(msg)

class ReadCloneFile(Target):
    ''' Read input clone file, return a "sample" obj
    '''
    def __init__(self, outdir, infile, parsefunc):
        Target.__init__(self)
        self.outdir = outdir
        self.infile = infile
        self.parsefunc = parsefunc

    def run(self):
        name = os.path.splitext(os.path.basename(self.outdir))[0]

        #split input file into smaller files if too big:
        global_dir = self.getGlobalTempDir()
        split_dir = os.path.join(global_dir, "samples_split", name)
        system("mkdir -p %s" % split_dir)
        lineperfile = 50000
        outbase = os.path.join(split_dir, name)
        split_file(self.infile, lineperfile, outbase)

        for f in os.listdir(split_dir):
            infile = os.path.join(split_dir, f)
            outfile = os.path.join(self.outdir, "%s.pickle" % f)
            self.addChildTarget(ReadCloneFileNoSplit(infile, outfile,
                                                               self.parsefunc))
        self.setFollowOnTarget(libcommon.CleanupDir(split_dir))

class ReadCloneFiles(Target):
    '''Set up children jobs to read each input file
    '''
    def __init__(self, outdir, indir, format="mitcr", ext=None):
        Target.__init__(self)
        self.outdir = outdir
        self.indir = indir
        self.format = format
        self.ext = ext

    def run(self):
        format2parsefunc = {"mitcr": mitcr.mitcr_parseline,
                            "adaptive": adaptive.adaptive_parseline,
                            "sequenta": sequenta.sequenta_parseline,
                            "aimseqtk": libclone.clone_parseline}
        if self.format not in format2parsefunc:
            raise UnrecognizedFormatError("Format %s is not recognized. Please\
                    choose one of these: %s\n" % 
                    (self.format, ",".join(format2parsefunc.keys())))

        parsefunc = format2parsefunc[self.format]

        for file in os.listdir(self.indir):
            # Check for the correct file extension if ext is specified
            basename = file
            if self.ext is not None:
                basename, extension = os.path.splitext(file)
                if not extension.endswith(self.ext):
                    continue
            infile = os.path.join(self.indir, file)
            samdir = os.path.join(self.outdir, "%s" % basename)
            system("mkdir -p %s" % samdir)
            self.addChildTarget(ReadCloneFile(samdir, infile, parsefunc))


