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
import numpy as np

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
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    group2line = {}
    lines = []
    linenames = []
    for name, obj in name2obj.iteritems():
        xdata, ydata = ld_get_sample_data(obj, attr)
        line, = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)
        if obj.group not in group2line:
            group2line[obj.group] = line

    if len(linenames) > 10:
        linenames = sorted(group2line.keys())
        lines = [group2line[group] for group in linenames]
    drawcommon.set_grid(axes)
    drawcommon.set_legend(axes, lines, linenames)
    drawcommon.edit_spine(axes)
    drawcommon.set_labels(axes, xlabel='Length', ylabel='%% of total %s' % attr)
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi) 

def draw_lendist_avr(name2obj, attr, outfile, outfmt='pdf', dpi=300):
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    g2x2y = {}
    g2numsam = {}
    g2color = {}
    minx = 30
    maxx = 0
    for name, obj in name2obj.iteritems():
        xdata, ydata = ld_get_sample_data(obj, attr)
        g = obj.group
        if g not in g2numsam:
            g2numsam[g] = 1
            g2color[g] = obj.color
        else:
            g2numsam[g] += 1
        minx = min(minx, min(xdata))
        maxx = max(maxx, max(xdata))
        for i, x in enumerate(xdata):
            y = ydata[i]
            if g not in g2x2y:
                g2x2y[g] = {x: [y]}
            elif x not in g2x2y[g]:
                g2x2y[g][x] = [y]
            else:
                g2x2y[g][x].append(y)

    lines = []
    linenames = []
    xdata = range(minx, maxx + 1)
    for g in sorted(g2x2y.keys()):
        numsam = g2numsam[g]
        mean_ydata = []
        std_ydata = []
        for x in xdata:
            if x not in g2x2y[g]:
                mean_ydata.append(0)
                std_ydata.append(0)
            else:
                ylist = g2x2y[g][x] + [0] * (numsam - len(g2x2y[g][x]))
                mean_ydata.append(np.mean(ylist))
                std_ydata.append(np.std(ylist))

        line, = axes.plot(xdata, mean_ydata, color=g2color[g], marker='o',
                         markeredgecolor=g2color[g], linestyle='-')
        axes.errorbar(xdata, mean_ydata, yerr=std_ydata, color=g2color[g],
                      linestyle='None', marker='None')
        lines.append(line)
        linenames.append(g)

    axes.set_ylim(bottom=-0.005)
    drawcommon.set_legend(axes, lines, linenames)
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_labels(axes, xlabel='Length', ylabel='%% of total %s' % attr)
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi) 

