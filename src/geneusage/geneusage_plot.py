#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw geneusage
'''

import os
import sys

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pyplot
from jobTree.scriptTree.target import Target

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

def draw_gene_usage(name2obj, attr, type, outbase, genes, opts):
    '''xaxis: genes
    yaxis: relative usage
    one line/ sample
    '''
    w = 10.0
    h = 8.0
    fig, pdf = drawcommon.init_image(w, h, opts.plotformat, outbase, opts.dpi)
    axes = drawcommon.set_axes(fig)

    genes = sorted(genes)
    xdata = range(len(genes))
    group2line = {}
    lines = []
    linenames = []
    for name, obj in name2obj.iteritems():
        ydata = gu_get_sample_data(obj, attr, type, genes)
        line = axes.plot(xdata, ydata, color=obj.color, marker=obj.marker,
                         markeredgecolor=obj.color, linestyle='-')
        lines.append(line)
        linenames.append(obj.name)
        if obj.group not in group2line:
            group2line[obj.group] = line
    
    if len(linenames) > 10:
        linenames = sorted(group2line.keys())
        lines = [group2line[group] for group in linenames]
    legend = axes.legend(lines, linenames, numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False

    drawcommon.edit_spine(axes)
    axes.xaxis.set_ticklabels(genes)
    axes.set_xlabel("Gene", size='x-large', weight='bold')
    axes.set_ylabel("% of total %s" % attr, size='x-large', weight='bold')
    drawcommon.write_image(fig, pdf, outformat, outbase, dpi) 
    
class GeneUsagePlot(Target):
    '''
    '''
    def __init__(self, name2obj, attr, type, genes, outbase, opts):
        Target.__init__(self)
        self.name2obj = name2obj
        self.attr = attr
        self.type = type
        self.genes = genes
        self.outbase = outbase
        self.opts = opts

    def run(self):
        if len(self.type) == 1:
            draw_gene_usage(self.name2obj, self.attr, self.type, self.genes,
                            self.outbase, self.opts)
        #else:  # draw matrix
        #    draw_combination_usage(self.name2obj, self.attr, self.type,
        #                           self.outbase, self.opts)
















