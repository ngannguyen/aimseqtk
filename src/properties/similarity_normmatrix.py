#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Compare sample's similarity
Input: a square matrix (rows=samples, cols=samples, cells=similarity)
       clonesize info (sample --> number of clones)
Output:
       normalized matrix
'''

import os
import sys
import re
import math
import numbers
import cPickle as pickle
import gzip
import numpy as np

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.common as lcommon

def read_clonesize(file):
    samples = []
    s2clones = {}

    f = open(file, 'r')
    f.readline()
    for line in f:
        items = line.strip().split('\t')
        sample = items[0]
        clones = int(items[1])
        if not re.search("Avr", sample):
            samples.append(sample)
        s2clones[sample] = clones
    f.close()
    return samples, s2clones

def read_matrix(file):
    s2s2c = {}
    f = open(file, 'r')
    colnames = f.readline().strip().split('\t')
    for line in f:
        items = line.strip().split('\t')
        rowname = items[0]
        s2s2c[rowname] = {}
        for i in xrange(1, len(items)):
            s2s2c[rowname][colnames[i]] = float(items[i])
    f.close()
    return s2s2c

def normalize_matrix(r2c2count, names, name2total):
    p_mean = 3.4 * 10**(-10)  # from Murugan 2012
    rows = []
    minval = float('inf')
    for rowname in names:
        row = []
        c2count = r2c2count[rowname]
        total1 = name2total[rowname]
        assert total1 > 0
        for colname in names:
            total2 = name2total[colname]
            assert total2 > 0
            newcount = c2count[colname]
            #newcount = c2count[colname] / (total1 * total2 * p_mean)
            #newcount = math.log10(c2count[colname] / (total1 * total2 * p_mean))
            minval = min(minval, newcount)
            row.append(newcount)
        rows.append(row)
    # adjust artifical values of self-overlap (so heatmap is visible)
    #minval = minval * 0.9
    for i, row in enumerate(rows):
        row[i] = minval
    return rows

def print_matrix(names, rows, file):
    f = open(file, 'w')
    f.write("Samples\t%s\n" % ("\t".join(names)))
    for i, name in enumerate(names):
        row = ["%f" % c for c in rows[i]]
        f.write("%s\t%s\n" % (name, "\t".join(row)))
    f.close()

def main():
    matrixfile = sys.argv[1]
    clonefile = sys.argv[2]
    outfile = sys.argv[3]
    names, name2total = read_clonesize(clonefile)
    r2c2count = read_matrix(matrixfile)
    matrix = normalize_matrix(r2c2count, names, name2total)
    print_matrix(names, matrix, outfile)

if __name__ == '__main__':
    main()


