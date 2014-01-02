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
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, opts.plotformat, outbase, opts.dpi)
    axes = drawcommon.set_axes(fig)

    xdata = range(len(groups))
    for row in rows:
        axes.plot(xdata, row, linestyle='-')
    
    drawcommon.edit_spine(axes)
    axes.xaxis.set_ticklabels(groups)
    axes.set_xlabel("Groups", size='x-large', weight='bold')
    axes.set_ylable("Clone size", size='x-large', weight='bold')
    axes.set_title("Clone %s" % clone, size='xx-large', weight='bold')
    drawcommon.write_image(fig, pdf, outformat, outbase, dpi)
