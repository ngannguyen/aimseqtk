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

def get_proper_hack(name):
    newname = name
    long2short = {"central-memory": "cMem", "naive": "N", "effector": "E",
                  "stemcell-memory": "scMem", "baseline": "BL",
                  "postchemo": "PC", "week6": "W6", "week14": "W14",
                  "week20": "W20"}
    for long, short in long2short.iteritems():
        newname = newname.replace(long, short)
    return newname

def draw_track_clone_hack(clone, rows, groups, outbase, opts):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    #sorted_groups = [
    #                 "cd4-naive-baseline", "cd4-naive-postchemo", "cd4-naive-week6", "cd4-naive-week14", "cd4-naive-week20",
    #                 "cd4-effector-baseline", "cd4-effector-postchemo", "cd4-effector-week6", "cd4-effector-week14", "cd4-effector-week20",
    #                 "cd4-central-memory-baseline", "cd4-central-memory-postchemo", "cd4-central-memory-week6", "cd4-central-memory-week14", "cd4-central-memory-week20",
    #                 "cd4-stemcell-memory-baseline", "cd4-stemcell-memory-postchemo", "cd4-stemcell-memory-week6", "cd4-stemcell-memory-week14", "cd4-stemcell-memory-week20",
    #                 "cd8-naive-baseline", "cd8-naive-postchemo", "cd8-naive-week6", "cd8-naive-week14", "cd8-naive-week20",
    #                 "cd8-effector-baseline", "cd8-effector-postchemo", "cd8-effector-week6", "cd8-effector-week14", "cd8-effector-week20",
    #                 "cd8-central-memory-baseline", "cd8-central-memory-postchemo", "cd8-central-memory-week6", "cd8-central-memory-week14", "cd8-central-memory-week20",
    #                 "cd8-stemcell-memory-baseline", "cd8-stemcell-memory-postchemo", "cd8-stemcell-memory-week6", "cd8-stemcell-memory-week14", "cd8-stemcell-memory-week20"
    #                 ]
    
    #sorted_groups = [
    #                 "cd4-N-BL", "cd4-N-PC", "cd4-N-W6", "cd4-N-W14", "cd4-N-W20",
    #                 "cd4-E-BL", "cd4-E-PC", "cd4-E-W6", "cd4-E-W14", "cd4-E-W20",
    #                 "cd4-cMem-BL", "cd4-cMem-PC", "cd4-cMem-W6", "cd4-cMem-W14", "cd4-cMem-W20",
    #                 "cd4-scMem-BL", "cd4-scMem-PC", "cd4-scMem-W6", "cd4-scMem-W14", "cd4-scMem-W20",
    #                 "cd8-N-BL", "cd8-N-PC", "cd8-N-W6", "cd8-N-W14", "cd8-N-W20",
    #                 "cd8-E-BL", "cd8-E-PC", "cd8-E-W6", "cd8-E-W14", "cd8-E-W20",
    #                 "cd8-cMem-BL", "cd8-cMem-PC", "cd8-cMem-W6", "cd8-cMem-W14", "cd8-cMem-W20",
    #                 "cd8-scMem-BL", "cd8-scMem-PC", "cd8-scMem-W6", "cd8-scMem-W14", "cd8-scMem-W20"
    #                 ]
    xdata = range(len(groups))
    #group2index = {}
    #for i, group in enumerate(groups):
    #    group2index[group] = i
    
    for row in rows:
        axes.plot(xdata, row, linestyle='-')
        #ydata = []
        #for group in sorted_groups:
        #    index = group2index[group]
        #    ydata.append(row[index])
        #axes.plot(xdata, ydata, linestyle='-')
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    #xlabels = []
    #for g in sorted_groups:
    #    xlabels.append(get_proper_hack(g))
    #drawcommon.set_xticks(axes, xdata, xlabels)
    #drawcommon.set_xticks(axes, xdata, sorted_groups)
    drawcommon.set_xticks(axes, xdata, groups)
    drawcommon.adjust_xticklabels(axes, size='xx-small', rotation=75)
    drawcommon.adjust_yticklabels(axes, size='xx-small')
    drawcommon.set_labels(axes, title="Clone %s" % clone, ylabel="Clone size")
    #drawcommon.set_labels(axes, "Clone %s" % clone, "Groups", "Clone size")
    
    drawcommon.write_image(fig, pdf, opts.plotformat, outbase, opts.dpi)

def draw_track_clone(clone, rows, groups, outbase, opts):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    xdata = range(len(groups))
    for row in rows:
        axes.plot(xdata, row, linestyle='-')
    
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, groups)
    drawcommon.adjust_xticklabels(axes, size='xx-small', rotation=75)
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

    
    
