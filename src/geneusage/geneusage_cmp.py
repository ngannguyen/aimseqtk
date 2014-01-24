#!/usr/bin/env python

import sys
import re
import os
import numpy as np

import aimseqtk.lib.drawcommon as drawcommon


def read_geneusage(file):
    g2n2genes = {}
    f = open(file, 'r')
    firstline = f.readline()
    items = firstline.strip().split('\t')
    if len(items) < 2:
        return
    genes = items[1:]
    n2genes = {}
    for line in f:
        items = line.strip().split('\t')
        n = items[0]
        if re.search('Avr', n):
            g2n2genes[n.split('_')[0]] = n2genes
            n2genes = {}
        else:
            n2genes[n] = [float(item) for item in items[1:]]
    f.close()
    return g2n2genes, genes

def get_median(vecs):
    if not vecs:
        return []
    med_vec = []
    l = len(vecs[0])
    for i in xrange(l):
        vec = [v[i] for v in vecs]
        med_vec.append(np.median(vec))
    return med_vec

def sample_diff(n2genes1, n2genes2):
    # compute the difference of each sample from condition 1 and condition2
    assert len(n2genes1) == len(n2genes2)
    n2diff = {}
    for n, genes1 in n2genes1.iteritems():
        genes2 = n2genes2[n]
        n2diff[n] = [g1 - genes2[i] for i, g1 in enumerate(genes1)]
    return n2diff

def sample_diff_med(n2genes1, n2genes2):
    # compute the difference of each sample from condition 1 and condition2
    assert len(n2genes1) == len(n2genes2)
    n2diff = {}
    med_vec = get_median(n2genes2.values())
    for n, genes1 in n2genes1.iteritems():
        n2diff[n] = [g1 - med_vec[i] for i, g1 in enumerate(genes1)]
    return n2diff

def union_genes(genes1, genes2):
    union_genes = genes1
    for g in genes2:
        if g not in union_genes:
            union_genes.append(g)
    return union_genes

def get_n2gene2freq(n2genes, genes):
    n2g2f = {}
    for n, ngenes in n2genes.iteritems():
        n2g2f[n] = {}
        for i, gene in enumerate(genes):
            n2g2f[n][gene] = ngenes[i]
    return n2g2f

def usage_union_genes(n2genes1, genes1, n2genes2, genes2):
    n2g2f_1 = get_n2gene2freq(n2genes1, genes1) 
    n2g2f_2 = get_n2gene2freq(n2genes2, genes2)
    genes = union_genes(genes1, genes2)
    n2genes1_union = {}
    n2genes2_union = {}
    for n, g2f_1 in n2g2f_1.iteritems():
        n2genes1_union[n] = []
        n2genes2_union[n] = []
        g2f_2 = n2g2f_2[n]
        for gene in genes:
            if gene in g2f_1:
                n2genes1_union[n].append(g2f_1[gene])
            else:
                n2genes1_union[n].append(0)
            if gene in g2f_2:
                n2genes2_union[n].append(g2f_2[gene])
            else:
                n2genes2_union[n].append(0)
    return n2genes1_union, n2genes2_union, genes

def diff_plot(n2diff, genes, outbase):
    axes, fig, pdf = drawcommon.get_axes(outfile=outbase)
    xdata = range(len(genes))
    data = []
    #for name, ydata in n2diff.iteritems():
        #axes.plot(xdata, ydata, linestyle='None', marker='.')
    for x in xdata:
        ydata = [n2diff[name][x] for name in n2diff]
        data.append(ydata)
    axes.boxplot(data)
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    xlabels = [g.lstrip("TRB") for g in genes]
    drawcommon.set_xticks(axes, [x + 1 for x in xdata], xlabels)
    axes.set_xlim(-0.5, len(xlabels) + 0.5)
    #axes.set_ylim(bottom=-0.005)
    drawcommon.adjust_ticklabels(axes, xrotation=75)
    drawcommon.set_labels(axes, xlabel="Gene", ylabel="Difference in gene usage")
    
    drawcommon.write_image(fig, pdf, 'pdf', outbase, 300) 

def pair_group_cmp(medvec, vecs, genes, outfile):
    axes, fig, pdf = drawcommon.get_axes(outfile=outfile)
    xdata = range(len(genes))
    data = []
    #for vec in vecs:
        #ydata = [v - medvec[i] for i, v in enumerate(vec)]
        #axes.plot(xdata, ydata, linestyle='None', marker='.')
    for x in xdata:
        ydata = [vec[x] - medvec[x] for vec in vecs]
        data.append(ydata)
    axes.boxplot(data)
    drawcommon.set_grid(axes)
    drawcommon.edit_spine(axes)
    xlabels = [g.lstrip("TRB") for g in genes]
    #drawcommon.set_xticks(axes, xdata, xlabels)
    drawcommon.set_xticks(axes, [x + 1 for x in xdata], xlabels)
    axes.set_xlim(-0.5, len(xlabels) + 0.5)
    #axes.set_ylim(bottom=-0.005)
    drawcommon.adjust_ticklabels(axes, xrotation=75)
    drawcommon.set_labels(axes, xlabel="Gene", ylabel="Difference in gene usage")
    
    drawcommon.write_image(fig, pdf, 'pdf', outfile, 300) 

def group_cmp(g2n2genes, genes, outbase):
    groups = g2n2genes.keys()
    for g1 in groups:
        n2genes1 = g2n2genes[g1]
        med1 = get_median(n2genes1.values())
        for g2 in groups:
            if g1 == g2:
                continue
            n2genes2 = g2n2genes[g2]
            outfile = "%s%s_%s" % (outbase, g1, g2)
            pair_group_cmp(med1, n2genes2.values(), genes, outfile)

def main():
    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    outdir = sys.argv[3]
    g2n2genes1, genes1 = read_geneusage(infile1) 
    g2n2genes2, genes2 = read_geneusage(infile2)
    for g, n2genes1 in g2n2genes1.iteritems():
        outfile = os.path.join(outdir, g)
        n2genes2 = g2n2genes2[g]
        if genes1 != genes2:
            n2genes1_union, n2genes2_union, genes_union = usage_union_genes(n2genes1, genes1, n2genes2, genes2)
        else:
            n2genes1_union = n2genes1
            n2genes2_union = n2genes2
            genes_union = genes1
        n2diff = sample_diff(n2genes1, n2genes2)
        #n2diff = sample_diff_med(n2genes1_union, n2genes2_union)
        diff_plot(n2diff, genes_union, outfile)

    # compare pair of groups:
    group_cmp(g2n2genes1, genes1, "%s/np_" % outdir)
    group_cmp(g2n2genes2, genes2, "%s/p_" % outdir)

if __name__ == '__main__':
    main()


