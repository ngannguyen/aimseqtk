#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw geneusage
'''

import os
import sys
import numpy as np

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot
from jobTree.scriptTree.target import Target

import aimseqtk.lib.common as libcommon
import aimseqtk.lib.drawcommon as drawcommon
import aimseqtk.lib.statcommon as statcommon


def gu_get_sample_data(stat, count_type, gene_type, genes):
    attr = "type2gene2%s" % count_type
    gene2count = stat[attr][gene_type]
    counts = []
    for gene in genes:
        if gene not in gene2count:
            counts.append(0.0)
        else:
            counts.append(gene2count[gene])
    return counts

def draw_pca(rownames, rows, outbase, name2obj, var1=None, var2=None):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase)
    lines = []
    linenames = []
    g2xdata = {}
    g2ydata = {}
    g2color = {}
    for i, (name, group) in enumerate(rownames):
        row = rows[i]
        x = row[0]
        y = row[1]
        if group not in g2xdata:
            g2xdata[group] = [x]
            g2ydata[group] = [y]
            g2color[group] = name2obj[name].color
        else:
            g2xdata[group].append(x)
            g2ydata[group].append(y)
    lines = []
    linenames = sorted(g2xdata.keys())
    for g in linenames:
        xdata = g2xdata[g]
        ydata = g2ydata[g]
        color = g2color[g]
        l = axes.plot(xdata, ydata, color=color, markersize=10.0, marker='.',
                      markeredgecolor=color, linestyle='none')
        lines.append(l)
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    if var1 and var2:
        drawcommon.set_labels(axes, xlabel="PC1, %f%%" % var1,
                              ylabel="PC2, %f%%" % var2)
    else:
        drawcommon.set_labels(axes, xlabel="PC1", ylabel="PC2")
    drawcommon.write_image(fig, pdf, outname=outbase) 

def draw_gene_usage(name2obj, attr, type, outbase, genes, opts=None):
    '''xaxis: genes
    yaxis: relative usage
    one line/ sample
    '''
    if opts:
        axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    else:
        axes, fig, pdf = drawcommon.get_axes(outfile=outbase)
    #genes = sorted(genes)
    genes = libcommon.sort_by_gene_number(genes)
    xdata = range(len(genes))
    group2line = {}
    lines = []
    linenames = []
    for name, obj in name2obj.iteritems():
        ydata = gu_get_sample_data(obj, attr, type, genes)
        line, = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)
        if obj.group not in group2line:
            group2line[obj.group] = line
    
    if len(linenames) > 10:
        linenames = sorted(group2line.keys())
        lines = [group2line[group] for group in linenames]
    
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.adjust_ticklabels(axes, xrotation=75)
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    xlabels = [g.lstrip("TRB") for g in genes]
    drawcommon.set_xticks(axes, xdata, xlabels)
    axes.set_xlim(-0.5, len(xlabels) + 0.5)
    axes.set_ylim(bottom=-0.005)
    drawcommon.adjust_ticklabels(axes, xrotation=75)
    drawcommon.set_labels(axes, xlabel="Gene", ylabel="%% of total %s" % attr)
   
    if opts:
        drawcommon.write_image(fig, pdf, opts.plotformat, outbase, opts.dpi) 
    else:
        drawcommon.write_image(fig, pdf, outname=outbase) 

def draw_gene_usage_avr(name2obj, attr, type, outbase, genes, opts, bar=False):
    '''xaxis: genes
    yaxis: relative usage
    one line/ sample
    '''
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase,
                                         outfmt=opts.plotformat, dpi=opts.dpi)
    #genes = sorted(genes)
    #genes = sorted(genes, key=lambda g: libcommon.get_gene_number(g))
    genes = libcommon.sort_by_gene_number(genes)
    xdata = range(len(genes))
    g2ydata = {}
    g2color = {}
    for name, obj in name2obj.iteritems():
        ydata = gu_get_sample_data(obj, attr, type, genes)
        g = obj.group
        if g not in g2ydata:
            g2ydata[g] = [[y] for y in ydata]
            g2color[g] = obj.color
        else:
            assert len(g2ydata[g]) == len(ydata)
            for i, y in enumerate(ydata):
                g2ydata[g][i].append(y)

    barwidth = (1.0 - 0.35) / len(g2ydata.keys())
    lines = []
    linenames = []
    #for g, ydata in g2ydata.iteritems():
    for i, g in enumerate(sorted(g2ydata.keys())):
        ydata = g2ydata[g]
        ydata = [ylist if ylist else [0.0] for ylist in ydata]
        mean_ydata = [np.mean(ylist) for ylist in ydata]
        std_ydata = [np.std(ylist) for ylist in ydata]
        if not bar:
            line, = axes.plot(xdata, mean_ydata, color=g2color[g],
                              linestyle='-', markeredgecolor=g2color[g],
                              marker='o')
            lines.append(line)
            axes.errorbar(xdata, mean_ydata, yerr=std_ydata, color=g2color[g],
                          linestyle="None", marker="None")
        else:
            group_xdata = [x + barwidth * i for x in xdata]
            line = axes.bar(group_xdata, mean_ydata, barwidth, yerr=std_ydata,
                            color=g2color[g], ecolor="#424242",
                            edgecolor=g2color[g])
            lines.append(line[0])
        linenames.append(g)
    
    drawcommon.set_legend(axes, lines, linenames) 
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    xlabels = [g.lstrip("TRB") for g in genes]
    if not bar:
        drawcommon.set_xticks(axes, xdata, xlabels)
    else: 
        xticks = [x + 0.325 for x in xdata]
        drawcommon.set_xticks(axes, xticks, xlabels)
    axes.set_xlim(-0.5, len(xlabels) + 0.5)
    axes.set_ylim(bottom=-0.005)
    drawcommon.adjust_ticklabels(axes, xrotation=75)
    drawcommon.set_labels(axes, xlabel="Gene", ylabel="%% of total %s" % attr)
    
    drawcommon.write_image(fig, pdf, opts.plotformat, outbase, opts.dpi) 
    
class GeneUsagePlot(Target):
    '''
    '''
    def __init__(self, name2obj, attr, type, genes, outbase, opts, avr=False):
        Target.__init__(self)
        self.name2obj = name2obj
        self.attr = attr
        self.type = type
        self.genes = genes
        self.outbase = outbase
        self.opts = opts
        self.avr = avr

    def run(self):
        if len(self.type) == 1:
            if self.avr or len(self.name2obj) > 20:
                draw_gene_usage_avr(self.name2obj, self.attr, self.type,
                                self.outbase, self.genes, self.opts, bar=True)
            else:
                draw_gene_usage(self.name2obj, self.attr, self.type,
                                self.outbase, self.genes, self.opts)
        #else:  # draw matrix
        #    draw_combination_usage(self.name2obj, self.attr, self.type,
        #                           self.outbase, self.opts)
















