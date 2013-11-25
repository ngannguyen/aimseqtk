#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw rarefaction plots
'''

import os
import sys

import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties

import aimseqtk.lib.drawcommon as drawcommon


def draw_rarefaction(name2size2sampling, name2sample, groups, index, outfile,
                     outformat='pdf', dpi=300):
    # xaxis: sampling size; yaxis: index value +- std
    width = 10.0
    height = 8.0
    fig, pdf = drawcommon.init_image(width, height, outformat, outfile, dpi)
    axes = drawcommon.set_axes(fig)
    
    lines = []
    linenames = []
    for name, size2sampling in name2size2sampling.iteritems():
        xdata = sorted(size2sampling.keys())
        ydata = []
        stddata = []
        for x in xdata:
            sampling = size2sampling[x]
            y = sampling[index]
            ydata.append(y)
            stdindex = "%s_std" % index
            if stdindex in sampling:
                std = sampling[stdindex]
                stddata.append(std)

        sample = name2sample[name]
        if stddata:
            axes.errorbar(xdata, ydata, yerr=stddata, color=sample.color,
                          markeredgecolor=color, fmt='.')
        line = axes.plot(xdata, ydata, color=sample.color, linestyle='-')
        if not groups:
            lines.append(line)
            linenames.append(sample.name)
        else:
            if sample.group not in linenames:
                lines.append(line)
                linenames.append(sample.group)
    
    axes.xaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    axes.yaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    drawcommon.edit_spine(axes)
    # Labeling:
    axes.set_title("%s Rarefaction Curve" % index.title(), size='xx-large', 
                   weight='bold')
    axes.set_xlabel("Sampling size (number of sequences)", size='x-large',
                   weight='bold')
    axes.set_ylabel(index.title(), size='x-large', weight='bold')
    
    legend = axes.legend(lines, linenames, numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False
    
    pyplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    pyplot.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outformat, outfile, dpi) 


