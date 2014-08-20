#!/usr/bin/env python

'''
Wed Jul  9 23:23:20 PDT 2014
Input: diversity.txt
Output: time-series plot
'''

import os
import sys
import re

from sonLib.bioio import system

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.drawcommon as drawcommon

def read_in_file(file, index="shannon"):
    sam2time2val = {}
    f = open(file, 'r')
    header_items = f.readline().split('\t')
    col = 0
    for i, header in enumerate(header_items):
        if header == 'shannon':
            col = i
            break
    if col == 0:
        sys.stderr.write("Wrong format\n")
        sys.exit(1)

    for line in f:
        items = line.rstrip().split('\t')
        names = items[0].split('-')
        sam = names[0]
        time = names[1]
        v = float(items[col])
        if sam not in sam2time2val:
            sam2time2val[sam] = {time: v}
        else:
            sam2time2val[sam][time] = v
    f.close()
    return sam2time2val

def draw_plot(sam2time2val, outfile, timepoints=None):
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile)
    
    group2color = {'IL7_baseline': "#c6dbef", 'IL7_postchemo': "#9ecae1",
                   'IL7_week6': "#6baed6", 'IL7_week14': "#3182bd",
                   'IL7_week20': "#08519c", 'none_baseline': "#d9d9d9",
                   'none_postchemo': "#bdbdbd", 'none_week6': "#969696",
                   'none_week14': "#636363", 'none_week20': "#252525"}
    groups = ['IL7', 'none']
    group2sam = {'IL7': ['LaBo', 'FoCh', 'JaMa'], 'none': ['LiBr', 'BrBu', 'MeRi']}
    if timepoints is None:
        timepoints = ['baseline', 'postchemo', 'week6', 'week14', 'week20']
    
    barwidth = (1.0 - 0.25) / len(timepoints)
    
    sam_per_group = 3  #  number of sample per group
    xticks = []
    
    ymin = 20
    ymax = 0

    for i, group in enumerate(groups):
        xdata = [x + i * sam_per_group + 1 for x in xrange(sam_per_group)]
        samples = group2sam[group]
        num_times = len(timepoints)
        group_xticks = [x + float(num_times) / 2.0 * barwidth for x in xdata]
        xticks.extend(group_xticks)

        for j, timepoint in enumerate(timepoints):
            t_xdata = [x + j * barwidth for x in xdata]
            color = group2color["%s_%s" % (group, timepoint)]
            ydata = []
            for sam in samples:
                y = sam2time2val[sam][timepoint]
                ydata.append(y)
                ymin = min(ymin, y)
                ymax = max(ymax, y)
            line = axes.bar(t_xdata, ydata, barwidth, color=color, ecolor="#424242")
  
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    xlabels = ['LaBo_IL7', 'FoCh_IL7', 'JaMa_IL7', 'LiBr_none', 'BrBu_none', 'MeRi_none']
    drawcommon.set_xticks(axes, xticks, xlabels)
    #axes.set_xlim(xmin=8, xmax=min(maxx, 30))
    axes.set_ylim(ymin - 1, ymax + 0.5)

    drawcommon.set_labels(axes, xlabel='Samples', ylabel='Shannon Diversity Index')
    drawcommon.write_image(fig, pdf, outname=outfile)
  
def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    timepoints = None
    if len(sys.argv) > 3:
        timepoints = sys.argv[3].split(',')
    sam2time2val = read_in_file(infile)
    draw_plot(sam2time2val, outfile, timepoints)
     

if __name__ == '__main__':
    main()

