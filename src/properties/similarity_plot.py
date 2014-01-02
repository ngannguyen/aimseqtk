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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def draw_group_pairwise_similarity(group1, group2, vec11, vec12, vec22,
                                   outfile, outfmt='pdf', dpi=300):
    # xaxis: G1_G1, G1_G2, G2_G2; yaxis: similarity index
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, outfmt, outfile, dpi)
    axes = drawcommon.set_axes(fig)

    xlabels = [group1.title(), '%s_%s' % (group1.title(), group2.title()),
               group2.title()]
    data = [vec11, vec12, vec22]
    axes.boxplot(data)
    drawcommon.edit_spine(axes)
    axes.xaxis.set_ticklabels(xlabels)

    #axes.set_xlabel("Group", size='x-large', weight='bold')
    axes.set_ylabel(attr.title(), size='x-large', weight='bold')
    pytplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)


