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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

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
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, outfmt, outfile, dpi)
    axes = drawcommon.set_axes(fig)

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
        line = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)

    legend = axes.legend(lines, linenames, numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False
    drawcommon.edit_spine(axes)
    axes.xaxis.set_ticklabels(xlabels)
    
    labels = cs_get_attr_plot_labels(attr, numtop)
    axes.set_title(labels[0], size='xx-large', weight='bold')
    axes.set_xlabel(labels[1], size='x-large', weight='bold')
    axes.set_ylabel(labels[2], size='x-large', weight='bold')

    drawcommon.write_image(fig, pdf, outformat, outfile, dpi) 

