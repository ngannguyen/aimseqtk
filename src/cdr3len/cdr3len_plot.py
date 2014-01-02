#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw CDR3 len distributions
xaxis = aa length
yaxis = % total clones OR % total reads
'''

import os
import sys

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def ld_get_sample_data(obj, attr):
    attr = "len2" + attr
    len2freq = obj[attr]
    xdata = sorted(len2freq.keys())
    ydata = [len2freq[x] for x in xdata]
    return xdata, ydata

def draw_lendist(name2obj, attr, outfile, outfmt='pdf', dpi=300):
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, outfmt, outfile, dpi)
    axes = drawcommon.set_axes(fig)

    group2line = {}
    lines = []
    linenames = []
    for name, obj in name2obj.iteritems():
        xdata, ydata = ld_get_sample_data(obj, attr)
        line = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)
        if obj.group not in group2line:
            group2line[group] = line

    if len(linenames) > 10:
        linenames = sorted(group2line.keys())
        lines = [group2line[group] for group in linenames]
    legend = axes.legend(lines, linenames, numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False
    
    drawcommon.edit_spine(axes)
    axes.set_xlabel("Length", size='x-large', weight='bold')
    axes.set_ylabel("% of total %s" % attr, size='x-large', weight='bold')
    drawcommon.write_image(fig, pdf, outformat, outfile, dpi) 

