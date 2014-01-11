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

    
