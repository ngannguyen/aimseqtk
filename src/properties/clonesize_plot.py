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

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def cs_get_attr_plot_labels(attr, numtop=50):
    attr2labels = {'numclones': ('Distribution of Clones',
                                 'Clone size (% of total sequences)', 
                                 'Frequency (% of total clones)'),
                   'counts': ('Distribution of Sequences',
                              'Clone size (% of total sequences)', 
                              'Frequency (% of total sequences)'),
                   'topfreqs': ('Distribution of %d Largest Clones' % numtop,
                                'Number of top clones',
                                'Clone Rank'),
                   'numclones_cumul': ('Cumulative Distribution of Clones',
                                 'Clone size (% of total sequences)', 
                                 'Frequency (% of total clones)'),
                   'counts_cumul': ('Cumulative Distribution of Sequences',
                              'Clone size (% of total sequences)', 
                              'Frequency (% of total sequences)'),
                   'topfreqs_cumul':
                      ('Cumulative Distribution of %d Largest Clones' % numtop,
                                'Number of top clones',
                                'Frequency (% of total sequences)')
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

