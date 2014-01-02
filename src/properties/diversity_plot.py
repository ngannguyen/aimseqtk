#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw diversity plots
xaxis = groups
yaxis = diversity index values
'''

import os
import sys

import aimseqtk.lib.drawcommon as drawcommon

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot


def draw_diversity_plot(group2names, name2obj, attr, outfile, outfmt='pdf',
                                                                      dpi=300):
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, outfmt, outfile, dpi)
    axes = drawcommon.set_axes(fig)
    
    xlabels = sorted(group2names.keys())
    data = []
    for group in xlabels:
        names = group2names[group]
        vec = [name2obj[name][attr] for name in names]
        data.append(vec)
    axes.boxplot(data)
    drawcommon.edit_spine(axes)
    axes.xaxis.set_ticklabels(xlabels)
    # Labeling:
    #axes.set_title("%s" % attr.title(), size='xx-large', weight='bold')
    axes.set_xlabel("Group", size='x-large', weight='bold')
    axes.set_ylabel(attr.title(), size='x-large', weight='bold')
    pytplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)

    



