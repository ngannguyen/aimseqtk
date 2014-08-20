#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

import os
import sys
from math import log10

from sonLib.bioio import system

import aimseqtk.src.overlap.trackclone_plot as tcp


def draw_plots(infile, outdir):
    mingroupsam = 5
    minfreq = -4  # 0.1 %

    f = open(infile, 'r')
    groups = f.readline().rstrip('\n').split('\t')[1:]
    for line in f:
        items = line.rstrip('\n').split('\t')
        assert len(items) == len(groups) + 1
        clone = items[0]
        outfile = os.path.join(outdir, clone)
        data = []
        signi = False
        for item in items[1:]:
            groupdata = []
            if item:
                groupdata = [log10(float(v)) for v in item.split(',')]
            data.append(groupdata)
            numpass_minfreq = 0
            for d in groupdata:
                if d >= minfreq:
                    numpass_minfreq += 1
            if numpass_minfreq >= mingroupsam:
                signi = True
        if signi:
            tcp.draw_track_clone_no_matched2(clone, data, groups, outfile) 
    f.close()

def main():
    infile = sys.argv[1]
    outdir = sys.argv[2]
    draw_plots(infile, outdir)    

if __name__ == '__main__':
    main()

