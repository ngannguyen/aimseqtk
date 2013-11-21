#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.sample as libsample
import aimseqtk.src.input.mitcr as mitcr
import aimseqtk.src.input.adaptive as adaptive
import aimseqtk.src.input.sequenta as sequenta


class FormatError(Exception):
    pass

class UnrecognizedFormatError(Exception):
    pass


def read_clone_file(file, parsefunc=mitcr.mitcr_parseline):
    samplename = os.path.splitext(os.path.basename(file))[0]
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

    return libsample.Sample(samplename, clones)

def read_clone_files(indir, parsefunc=mitcr.mitcr_parseline, ext=None):
    # "func" is the converting function, different input formats
    name2sample = {}
    
    for file in os.listdir(indir):
        # Check for the correct file extension if ext is specified
        if ext is not None:
            basename, extension = os.path.splitext(file)
            if not extension.endswith(ext):
                continue

        sample = read_clone_file(os.path.join(indir, file), parsefunc)
        if sample.name in name2sample:
            sys.stderr.write("Warning: Repetitive sample %s\n" % sample.name)
            name2sample[sample.name].addclones(sample.clones)
        else:
            name2sample[sample.name] = sample
    
    return name2sample
        
def process_clone_inputs(indir, format="mitcr", ext=None):
    # format has to be one of the following values: mitcr, adaptive, sequenta
    # return name2sample (key = sample.name, val = Sample())
    if format == "mitcr":
        return read_clone_files(indir, ext=ext)
    elif format == "adaptive":
        return read_clone_files(indir, adaptive.adaptive_parseline, ext)
    elif format == "sequenta":
        return read_clone_files(indir, sequenta.sequenta_parseline, ext)
    else:
        types = ["mitcr", "sequenta", "adaptive"]
        raise UnrecognizedFormatError("Format %s is not recognized. Please \
                    choose one of these: %s\n" % (format, ",".join(types)))
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
        if group in group2samples:
            sys.stderr.write("Warning: group_info file %s has groups with the "
                             + "same name: %s -- merge them into 1 group.\n" % 
                             (file, group))
            group2samples[group].extend(samples)
        else:
            group2samples[group] = samples

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

    return group2samples









