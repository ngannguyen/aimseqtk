#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw clone frequencies over different groups
xaxis: groups (e.g different time points)
yaxis: clone size (frequency)
one line per sample
'''

import os
import sys
import math

import aimseqtk.lib.drawcommon as drawcommon

def draw_track_clone(clone, rows, groups, outbase, opts):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    xdata = range(len(groups))
    for row in rows:
        axes.plot(xdata, row, linestyle='-')
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, groups)
    drawcommon.set_labels(axes, "Clone %s" % clone, "Groups", "Clone size")
    
    drawcommon.write_image(fig, pdf, opts.plotformat, outbase, opts.dpi)

def draw_track_clone_no_matched(clone, data, groups, outbase, opts):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    xdata = range(1, len(groups) + 1)
    axes.boxplot(data)
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, groups)
    #drawcommon.set_labels(axes, "Clone %s" % clone, "Groups", "Clone size")
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, opts.plotformat, outbase, opts.dpi)

def draw_track_clone_no_matched2(clone, data, groups, outbase, opts=None):
    if opts is not None:
        axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    else:
        axes, fig, pdf = drawcommon.get_axes(outfile=outbase)
    xdata = range(1, len(groups) + 1)
    offset = 0.02
    numpoints = 25.0
    for i, x in enumerate(xdata):
        ydata = sorted(data[i])
        halfwidth = min(numpoints, len(ydata)) * offset / 2
        currxdata = []
        numrounds = int(math.ceil(len(ydata) / numpoints))
        for r in xrange(numrounds):
            numy = min(numpoints, len(ydata) - r * numpoints)
            for j in xrange(int(numy)):
                currxdata.append(x - halfwidth + j * offset)
        axes.plot(currxdata, ydata, ls='none', marker='o', markersize=5)

    #axes.boxplot(data)
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, groups)
    #drawcommon.set_labels(axes, "Clone %s" % clone, "Groups", "Clone size")
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axes.set_xlim(0, len(groups) + 1)

    drawcommon.write_image(fig, pdf, outname=outbase)

    
    
