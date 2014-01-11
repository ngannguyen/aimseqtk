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
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    xdata = range(1, len(group2names) + 1)
    xlabels = sorted(group2names.keys())
    data = []
    for group in xlabels:
        names = group2names[group]
        vec = [name2obj[name][attr] for name in names]
        data.append(vec)
    axes.boxplot(data)
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, xlabels)
    drawcommon.set_labels(axes, xlabel="Group", ylabel=attr.title())
    
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)

    



