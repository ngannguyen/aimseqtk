#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys

import aimseqtk.lib.common as libcommon


class FormatError(Exception):
    pass

def read_input_file(parsefunc, file):
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
                           line for appropriate column names" % file)

    return Sample(samplename, clones)

def convert_input_files(parsefunc, indir, ext=None):
    # "func" is the converting function, different input formats
    name2sample = {}
    
    for file in os.listdir(indir):
        # Check for the correct file extension if ext is specified
        if ext is not None:
            basename, extension = os.path.splitext(file)
            if not extension.endswith(ext):
                continue

        sample = read_input_file(parsefunc, os.path.join(indir, file))
        if sample.name in name2sample:
            sys.stderr.write("Warning: Repetitive sample %s\n" % sample.name)
            name2sample[sample.name].addclones(sample.clones)
        else:
            name2sample[sample.name] = sample
    
    return name2sample
        


