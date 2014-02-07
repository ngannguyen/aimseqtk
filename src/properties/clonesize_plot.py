#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw clonesize plots
xaxis = clonesize categories (in freqs) or clone rank
yaxis = % total clones or % total seqs
cumulative or discrete
'''

import os
import sys
import numpy as np
from math import log10

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def cs_get_attr_plot_labels(attr, numtop=50):
    attr2labels = {'numclones': ('Distribution of Clones',
                                 'Clone size (proportion of total sequences)', 
                                 'Frequency (% of total clones), log10'),
                   'counts': ('Distribution of Sequences',
                              'Clone size (proportion of total sequences)', 
                              'Frequency (% of total sequences)'),
                   'topfreqs': ('Distribution of %d Largest Clones' % numtop,
                                'Clone rank',
                                'Clone frequency (% of total sequences)'),
                   'numclones_cumul': ('Cumulative Distribution of Clones',
                                 'Clone size (proportion of total sequences)', 
                                 'Frequency (% of total clones), log10'),
                   'counts_cumul': ('Cumulative Distribution of Sequences',
                              'Clone size (proportion of total sequences)', 
                              'Frequency (% of total sequences)'),
                   'topfreqs_cumul':
                      ('Cumulative Distribution of %d Largest Clones' % numtop,
                                'Number of top clones',
                                'Clone frequency (% of total sequences)')
                  }
    return attr2labels[attr]

def draw_clonesize_dist(name2obj, attr, outfile, outfmt='pdf', dpi=300):
    # name2stat: key = sample name, val = CloneSizeStat
    # attr = numclones/counts/numclones_cumul/counts_cumul/topfreqs/
    #         topfreqs_cumul
    assert len(name2obj) > 0
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    obj0 = name2obj.values()[0]
    xlabels = [str(x) for x in obj0.freqs]
    if attr == 'topfreqs' or attr == 'topfreqs_cumul':
        xlabels = [str(x + 1) for x in xrange(len(obj0.topfreqs))] 
    xdata = range(0, len(xlabels))
    numtop = len(xdata)
    if attr == 'numclones' or attr == 'numclones_attr':
        axes.set_yscale('log')
    linenames = []
    lines = []
    for name, obj in name2obj.iteritems():
        ydata = obj[attr]
        if len(xdata) != len(ydata):  # HACK
            continue
        line, = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)

    drawcommon.set_grid(axes)
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, xlabels)
    
    labels = cs_get_attr_plot_labels(attr, numtop)
    drawcommon.set_labels(axes, labels[0], labels[1], labels[2])

    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi) 

def get_logstd(ydata):
    std_ydata = []
    for ylist in ydata:
        log_ylist = []
        for y in ylist:
            if y > 0:
                log_ylist.append(log10(y))
        if log_ylist:
            std_ydata.append(np.std(log_ylist))
        else:
            std_ydata.append(0.0)
    return std_ydata
                
def draw_clonesize_dist_avr(name2obj, attr, outfile, outfmt='pdf', dpi=300):
    # name2stat: key = sample name, val = CloneSizeStat
    # attr = numclones/counts/numclones_cumul/counts_cumul/topfreqs/
    #         topfreqs_cumul
    assert len(name2obj) > 0
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    obj0 = name2obj.values()[0]
    xlabels = [str(x) for x in obj0.freqs]
    if attr == 'topfreqs' or attr == 'topfreqs_cumul':
        xlabels = [str(x + 1) for x in xrange(len(obj0.topfreqs))] 
    xdata = range(0, len(xlabels))
    numtop = len(xdata)
    #if attr == 'numclones':   # or attr == 'numclones_cumul':
    #    axes.set_yscale('log')
    
    g2ydata = {}
    g2color = {}
    for name, obj in name2obj.iteritems():
        ydata = obj[attr]
        g = obj.group
        if g not in g2ydata:
            g2ydata[g] = [[y] for y in ydata]
            g2color[g] = obj.color
        else:
            assert len(g2ydata[g]) == len(ydata)
            for i, y in enumerate(ydata):
                g2ydata[g][i].append(y)

    linenames = []
    lines = []
    #boxwidth = 0.05
    #offset = boxwidth + 0.01
    offset = 0.05
    for i, g in enumerate(sorted(g2ydata.keys())):
        g_xdata = [x + 0.5 + offset * (i - 1) for x in xdata]
        g_xdata[0] += 0.35
        ydata = g2ydata[g]
        #axes.boxplot(ydata, positions=g_xdata, widths=boxwidth)
        mean_ydata = [np.mean(ylist) for ylist in ydata]
        std_ydata = [np.std(ylist) for ylist in ydata]
        if attr == "numclones_cumul" or attr == "numclones":
            mean_ydata = [log10(y) if y > 0 else 2.1 for y in mean_ydata]
            std_ydata = get_logstd(ydata)
        line, = axes.plot(g_xdata, mean_ydata, color=g2color[g], linestyle='-',
                          markeredgecolor=g2color[g], marker='o', lw=2)
        lines.append(line)
        linenames.append(g)
        axes.errorbar(g_xdata, mean_ydata, yerr=std_ydata, color=g2color[g],
                      linestyle="None", marker="None")

    drawcommon.set_grid(axes)
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, xlabels)
    axes.set_xlim(-0.5, len(xlabels) + 0.5)
    if attr != "numclones_cumul" and attr != 'numclones':
        axes.set_ylim(bottom=-0.005)
    
    labels = cs_get_attr_plot_labels(attr, numtop)
    drawcommon.set_labels(axes, labels[0], labels[1], labels[2])

    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi) 

