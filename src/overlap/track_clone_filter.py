#!/usr/bin/env python

'''
Filter for clones that:
1/ Present at large frequency at baseline, disappear after chemo, but come back later
2/ Move from one population (e.g. effector) at an earlier timepoint to another population (e.g. memory) at a later timepoint
3/ Just appear after chemo
4/ Present at large frequency at baseline and after chemo
'''

import os
import sys

from sonLib.bioio import system



def read_infile(file):
    clone2sam2freq = {}
    f = open(file, 'r')
    header = f.readline()
    header_items = header.lstrip('#').rstrip('\n').split("\t")
    i2field = {}
    for i, field in enumerate(header_items):
        i2field[i] = field

    for line in f:
        items = line.rstrip('\n').split('\t')
        clone = items[0]
        sam2freq = {}
        for i in xrange(2, len(items)):
            freq = float(items[i])
            if freq > 0:
                sam = i2field[i]
                sam2freq[sam] = freq
        #if len(sam2freq) >= 2:  # clone is present in at least 2 samples
        if len(sam2freq) >= 10:  # clone is present in at least 8 samples
            clone2sam2freq[clone] = sam2freq

    return clone2sam2freq

def time_cmp(time1, time2):
    timepoints = {'BL': 0, 'PC': 1, 'W06': 2, 'W14': 3, 'W20': 4}
    return cmp(timepoints[time1], timepoints[time2])

def trackclone_filter(clone2sam2freq):
    '''
    1/ Present at large frequency at baseline, disappear after chemo, but come
       back later
    2/ Move from one population (e.g. effector) at an earlier timepoint to
       another population (e.g. memory) at a later timepoint
    3/ Just appear after chemo
    '''
    cat2clones = {1: [], 2: [], 3: []}
    for clone, sam2freq in clone2sam2freq.iteritems():
        cd2time2type = {}
        cd2baseline = {'cd4': False, 'cd8': False}

        for sam, freq in sam2freq.iteritems():
            items = sam.split('-')
            cd = items[0]
            if cd not in cd2time2type:
                cd2time2type[cd] = {}
            celltype = items[1]
            timepoint = items[2]
            if timepoint == 'BL' and freq >= 0.005:
                cd2baseline[cd] = True

            if timepoint not in cd2time2type[cd]:
                cd2time2type[cd][timepoint] = [celltype]
            else:
                cd2time2type[cd][timepoint].append(celltype)

        # find out which category this clone belongs to
        cats = []
        for cd, time2type in cd2time2type.iteritems():
            if len(time2type) == 1:
                continue
            times = sorted(time2type.keys(), cmp=time_cmp)
            if cd2baseline[cd] and 'BL' in times and 'PC' not in times:
                if 1 not in cats:
                    cats.append(1)
            for i in xrange(len(times) -1):
                time1 = times[i]
                types1 = time2type[time1]
                for j in xrange(i + 1, len(times)):
                    time2 = times[j]
                    types2 = time2type[time2]
                    if (("N" in types1 or "E" in types1) and
                        ("cMem" in types2 and "scMem" in types2)):
                       if 2 not in cats:
                           cats.append(2)
            if 'BL' not in times and 3 not in cats:
                cats.append(3)
        for c in cats:
            cat2clones[c].append(clone)

    return cat2clones

def print_results(cat2clones, outfile):
    f = open(outfile, 'w')
    f.write("#Summary: ")
    for c, clones in cat2clones.iteritems():
        f.write("%d: %d; " % (c, len(clones)))
    f.write("\n")
    for c, clones in cat2clones.iteritems():
        f.write("#%d\n" % c)
        for clone in clones:
            f.write("%s\n" % clone)
    f.close()

def separate_plots(cat2clones):
    outdir = "filter_10"
    system("mkdir -p %s" % outdir)
    for cat, clones in cat2clones.iteritems():
        catdir = os.path.join(outdir, "%d" % cat)
        system("mkdir -p %s" % catdir)
        for clone in clones:
            clonepath = os.path.abspath(os.path.join("plots", "%s.pdf" % clone))
            outpath = os.path.join(catdir, "%s.pdf" % clone)
            system("ln -s %s %s" % (clonepath, outpath))

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]

    clone2sam2freq = read_infile(infile)
    cat2clones = trackclone_filter(clone2sam2freq)
    print_results(cat2clones, outfile)
    separate_plots(cat2clones)

if __name__ == '__main__':
    main()

