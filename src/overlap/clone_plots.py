
'''
Draw geneusage and lendist 
Input: 
    lists of clones with format: V_cdr3seq_J, one file per group
'''

import sys
import os

from aimseqtk.lib.drawcommon import get_colors_medium as getcolors 
from aimseqtk.lib.common import union_lists
from aimseqtk.lib.statcommon import SampleStat
from aimseqtk.src.geneusage.geneusage_plot import draw_gene_usage
from aimseqtk.src.cdr3len.cdr3len_plot import draw_lendist

class CloneStat(SampleStat):
    def __init__(self, v2freq, j2freq, len2freq, color, marker, name):
        SampleStat.__init__(self)
        self.type2gene2clones = {}
        self.type2gene2clones['v'] = v2freq
        self.type2gene2clones['j'] = j2freq
        self.len2clones = len2freq
        self.color = color
        self.marker = marker
        self.name = name

def updatecount(mydict, k):
    if k in mydict:
        mydict[k] += 1
    else:
        mydict[k] = 1

def normalize(mydict):
    total = sum(mydict.values())
    for k, v in mydict.iteritems():
        mydict[k] = v*100.0/total

def readfile(file):
    v2freq = {}
    j2freq = {}
    len2freq = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split('_')
        v = items[0]
        seq = items[1]
        j = items[2]
        updatecount(v2freq, v)
        updatecount(j2freq, j)
        updatecount(len2freq, len(seq))
    normalize(v2freq)
    normalize(j2freq)
    normalize(len2freq)
    f.close()
    return v2freq, j2freq, len2freq

def readfiles(indir):
    name2obj = {}
    colors = getcolors()
    i = 0
    for file in os.listdir(indir):
        filepath = os.path.join(indir, file)
        v2f, j2f, l2f = readfile(filepath)
        name2obj[file] = CloneStat(v2f, j2f, l2f, colors[i], "o", file)
        i += 1
    return name2obj

def draw_plots(name2obj, outdir):
    vfile = os.path.join(outdir, "v")
    vgenes = union_lists([obj.type2gene2clones['v'].keys() for obj in name2obj.values()])
    draw_gene_usage(name2obj, "clones", "v", vfile, vgenes)

    jfile = os.path.join(outdir, "j")
    jgenes = union_lists([obj.type2gene2clones['j'].keys() for obj in name2obj.values()])
    draw_gene_usage(name2obj, "clones", "j", jfile, jgenes)

    lfile = os.path.join(outdir, "lendist")
    draw_lendist(name2obj, "clones", lfile)

def main():
    indir = sys.argv[1]
    outdir = sys.argv[2]
    name2obj = readfiles(indir)
    draw_plots(name2obj, outdir)
    
if __name__ == '__main__':
    main()

