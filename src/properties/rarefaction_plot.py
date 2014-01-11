#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw rarefaction plots
'''

import os
import sys

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot

import aimseqtk.lib.drawcommon as drawcommon


def draw_rarefaction(name2size2sampling, groups, index, outfile,
                     outformat='pdf', dpi=300):
    # xaxis: sampling size; yaxis: index value +- std
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outformat,
                                                                       dpi=dpi)

    lines = []
    linenames = []
    for name, size2sampling in name2size2sampling.iteritems():
        if not size2sampling:
            continue
        xdata = sorted(size2sampling.keys())
        ydata = []
        stddata = []
        for x in xdata:
            sampling = size2sampling[x]
            y = sampling[index]
            ydata.append(y)
            stdindex = "%s_std" % index
            if stdindex in sampling.getitems():
                std = sampling[stdindex]
                stddata.append(std)
        sampling0 = size2sampling[xdata[0]]
        color = sampling0.color
        marker = sampling0.marker
        if stddata:
            axes.errorbar(xdata, ydata, yerr=stddata, color=color,
                          markeredgecolor=color, fmt='.')
        line, = axes.plot(xdata, ydata, color=color, mec=color,
                         marker=marker, linestyle='-')
        if not groups:
            lines.append(line)
            linenames.append(sampling0.name)
        else:
            if sampling0.group not in linenames:
                lines.append(line)
                linenames.append(sampling0.group)
    
    drawcommon.set_grid(axes)
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.edit_spine(axes)
    # Labeling:
    title = "%s Rarefaction Curve" % index.title()
    xlabel = "Sampling size (number of sequences)"
    ylabel = index.title()
    drawcommon.set_labels(axes, title, xlabel, ylabel)
    
    axes.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outformat, outfile, dpi) 


