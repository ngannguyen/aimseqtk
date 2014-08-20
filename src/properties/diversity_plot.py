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
    drawcommon.adjust_ticklabels(axes, xrotation=30)
    drawcommon.set_labels(axes, xlabel="Group", ylabel=attr.title())
    
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)

def draw_diversity_plot_hacktimeseries(group2names, name2obj, attr, outfile,
                                       outfmt='pdf', dpi=300):
    '''Instead of box plot for each group, keep each sample separately
    '''
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile, outfmt=outfmt,
                                                                      dpi=dpi)
    xlabels = ['BL', 'PC', 'W6', 'W14', 'W20']
    xdata = range(1, len(xlabels) + 1)
    
    cat2group2names = {'None': {}, 'IL7': {}}
    for catgroup, names in group2names.iteritems():
        items = catgroup.split('-')
        cat = items[0]
        group = items[1]
        cat2group2names[cat][group] = names
    
    #HACK
    #sam2color = {'LiBr': "#6baed6", 'BrBu': "#3182bd", 'MeRi': "#08519c",  # blue
    #             'LaBo': "#74c476", 'FoCh': "#31a354", 'JaMa': "#006d2c"} 
    sam2color = {'LiBr': "#525252", 'BrBu': "#969696", 'MeRi': "#cccccc",  # gray
                 'LaBo': "#2171b5", 'FoCh': "#6baed6", 'JaMa': "#bdd7e7"}  # blue
    sam2data = {}
    sam2marker = {}
    name2cat = {}
    for cat, g2n in cat2group2names.iteritems():
        for i, group in enumerate(xlabels):  # each timepoint
            names = g2n[group]
            if i == 0:
                name2cat[names[0].split('-')[0]] = cat
            for name in names:
                sam = name.split('-')[0]
                y = name2obj[name][attr]
                if sam not in sam2data:
                    sam2data[sam] = [y]
                    #sam2color[sam] = name2obj[name].color
                    sam2marker[sam] = name2obj[name].marker
                else:
                    sam2data[sam].append(y)
    
    lines = []
    linenames = []
    for sam, ydata in sam2data.iteritems():
        l, = axes.plot(xdata, ydata, color=sam2color[sam],
                      marker=sam2marker[sam], linestyle='-',
                      markeredgecolor=sam2color[sam])
        if sam in name2cat:  #legend
            lines.append(l)
            linenames.append(name2cat[sam])
     
    drawcommon.set_grid(axes)
    drawcommon.set_legend(axes, lines, linenames)
    drawcommon.edit_spine(axes)
    drawcommon.set_xticks(axes, xdata, xlabels)
    drawcommon.adjust_ticklabels(axes, xrotation=30)
    axes.set_xlim(0.5, len(xdata) + 0.5)
    drawcommon.set_labels(axes, xlabel="Group", ylabel=attr.title())
    axes.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    drawcommon.write_image(fig, pdf, outfmt, outfile, dpi)
    



