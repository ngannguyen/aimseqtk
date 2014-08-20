#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Input: the similarity matrix and the group2samples file
Create heatmap without outliers. Create distribution of the similarities.
Report outliers
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
import aimseqtk.lib.drawcommon as dcommon
from aimseqtk.src.input.inputcommon import read_group_info
from aimseqtk.lib.statcommon import get_outliers, pair_vec
from aimseqtk.src.properties.similarity import table_group_pairwise_similarity
from aimseqtk.src.properties.similarity_plot import\
                                                draw_group_pairwise_similarity

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

def print_outliers(low_outliers, up_outliers, lowpoint, uppoint, outfile):
    f = open(outfile, 'w')
    f.write("#Name\tCount\n")
    f.write("#Upper bound is %s\n" % str(uppoint))
    for name, count in up_outliers.iteritems():
        f.write("%s\t%d\n" % (name, count))
    f.write("#Lower bound is %s\n" % str(lowpoint))
    for name, count in low_outliers.iteritems():
        f.write("%s\t%d\n" % (name, count))
    f.close()

def get_pair2count(r2c2count):
    pair2count = {}
    visited = {}
    for r in r2c2count.keys():
        visited[r] = []
    for r, c2count in r2c2count.iteritems():
        for c, count in c2count.iteritems():
            if r == c or (c in visited and r in visited[c]):
                continue
            visited[r].append(c)
            pair = "%s,%s" % (r, c)
            pair2count[pair] = count
    return pair2count

def get_dist(pair2count, g2s, outdir):
    low_outliers, up_outliers, lowpoint, uppoint = get_outliers(pair2count, 5)
    outliers_file = os.path.join(outdir, "outliers.txt")
    print_outliers(low_outliers, up_outliers, lowpoint, uppoint, outliers_file)
    
    dist_file = os.path.join("dist")
    #dist_plot(pair2count.values(), dist_file)
    return lowpoint, uppoint

def draw_heatmap(r2c2count, s2color, lowpoint, uppoint, outbase, names=None):
    if names is None:
        names = r2c2count.keys()
    rows = []
    for r in names:
        row = [r2c2count[r][c] for c in names]
        rows.append(row)
    dcommon.draw_heatmap(names, names, rows, outbase, minval=lowpoint,
                         maxval=uppoint, name2color=s2color, symm=True,
                         rcluster=False, ccluster=False)

def remove_outliers(vec, vmin, vmax):
    newvec = []
    for v in vec:
        if vmin <= v and v <= vmax:
            newvec.append(v)
    return newvec

def group_cmp(pair2count, g2s, ymin, ymax, outdir):
    groups = sorted(g2s.keys())
    if len(groups) < 2:
        return
    for i in xrange(0, len(groups) - 1):
        g1 = groups[i]
        names1 = g2s[g1]
        vec11 = pair_vec(names1, names1, pair2count, sep=',')
        vec11_ = remove_outliers(vec11, ymin, ymax)
        for j in xrange(i + 1, len(groups)):
            g2 = groups[j]
            names2 = g2s[g2]
            vec12 = pair_vec(names1, names2, pair2count, sep=',')
            vec22 = pair_vec(names2, names2, pair2count, sep=',')
            plotfile = os.path.join(outdir, "%s_%s" % (g1, g2))
            draw_group_pairwise_similarity(g1, g2, vec11, vec12,
                                         vec22, plotfile, ymin=ymin, ymax=ymax)
            tabfile = os.path.join(outdir, "%s_%s.txt" % (g1, g2))
            table_group_pairwise_similarity(g1, g2, vec11, vec12, vec22,
                                            tabfile, "ttest")
            # remove outliers:
            vec12_ = remove_outliers(vec12, ymin, ymax)
            vec22_ = remove_outliers(vec22, ymin, ymax)
            tabfile = os.path.join(outdir, "rmoutliers_%s_%s.txt" % (g1, g2))
            table_group_pairwise_similarity(g1, g2, vec11_, vec12_, vec22_,
                                            tabfile, "ttest")


def main():
    matrixfile = sys.argv[1]
    g2sfile = sys.argv[2]
    outdir = sys.argv[3]
    # read in matrix:
    r2c2count = read_matrix(matrixfile)
    # group info:
    groups, g2s, matched = read_group_info(g2sfile)

    pair2count = get_pair2count(r2c2count)
    # get distribution (draw dist and get outliers)
    lowpoint, uppoint = get_dist(pair2count, g2s, outdir)
    
    # re draw group comparison boxplot with lims:
    group_cmp(pair2count, g2s, lowpoint, uppoint, outdir)

    # draw heatmap:
    s2g = lcommon.get_val2key_1to1(g2s)
    s2color = dcommon.get_name2color_wtgroup(s2g.keys(), s2g, g2s)
    heatmap_file = os.path.join(outdir, "heatmap.pdf")
    names = []
    for g, gnames in g2s.iteritems():
        for n in gnames:
            if n in r2c2count:
                names.append(n)
    draw_heatmap(r2c2count, s2color, lowpoint, uppoint, heatmap_file, names) 

if __name__ == '__main__':
    main()


