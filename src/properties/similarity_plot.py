#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw similarity plots
1/ Heatmaps
2/ Group comparisons plot:
    xaxis: G1_G1, G1_G2, G2_G2
    yaxis: similarity indices of pair of samples
'''

import os
import sys

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def draw_group_pairwise_similarity(group1, group2, vec11, vec12, vec22,
                                   outfile, outfmt='pdf', dpi=300,
                                   xmax=None, xmin=None, ymax=None, ymin=None):
    # xaxis: G1_G1, G1_G2, G2_G2; yaxis: similarity index
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    xdata = range(1, 4)
    xlabels = [group1, '%s_%s' % (group1, group2), group2]
    data = [vec11, vec12, vec22]
    axes.boxplot(data)
    
    #axes.yaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, xlabels)

    #axes.set_xlabel("Group", size='x-large', weight='bold')
    #axes.set_ylabel(attr, size='x-large', weight='bold')
    if not ymax or ymax > 1000:
        axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    if xmin and xmax:
        axes.set_xlim(xmin=xmin, xmax=xmax)
    else:
        if xmin:
            axes.set_xlim(xmin=xmin)
        elif xmax:
            axes.set_xlim(xmax=xmax)
    if ymin and ymax:
        axes.set_ylim(ymin=ymin, ymax=ymax)
    else:
        if ymin:
            axes.set_ylim(ymin=ymin)
        elif ymax:
            axes.set_ylim(ymax=ymax)
    
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)


