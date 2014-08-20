#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''recombination prob. plots
'''

import os
import sys
import numbers
import cPickle as pickle
import gzip
import numpy as np

import aimseqtk.lib.common as lcommon
import aimseqtk.lib.drawcommon as drawcommon
from aimseqtk.src.recomb.recomb_model import RecombModel
import aimseqtk.src.recomb.recomb_cmp as rcmp
import aimseqtk.src.recomb.recomb_common as rcommon


def model_get_ydata(models, attr):
    assert attr in models[0].get_attrs_depth1()
    k2ydata = {}
    for k in models[0][attr]:
        k2ydata[k] = [m[attr][k] for m in models]
    return k2ydata

def model_get_ydata2(models, attr):
    assert attr in models[0].get_attrs_depth2()
    k_to_k2ydata = {}
    for k, k2v in models[0][attr].iteritems():
        k2ydata = {}
        for k2 in k2v:
            k2ydata[k2] = [m[attr][k][k2] for m in models]
        k_to_k2ydata[k] = k2ydata
    return k_to_k2ydata

def model_makeplots(indir, outdir):
    objs = lcommon.load_pickledir(indir)
    if not objs:
        return
    rcommon.union_models(objs)
    for attr in objs[0].get_attrs_depth1():
        outbase = os.path.join(outdir, attr)
        k2ydata = model_get_ydata(objs, attr)
        gene_attrs = ['v', 'd', 'j']
        if attr in gene_attrs:
            xlabels = lcommon.sort_by_gene_number(k2ydata.keys())
            rcmp.diff_plot(k2ydata, outbase, xlabels, attr.upper(),
                           ylabel='Frequency')
        else:
            rcmp.diff_plot(k2ydata, outbase, label=attr.upper(), xmin=-0.5,
                           xmax=30.5, ylabel='Frequency')
    for attr in objs[0].get_attrs_depth2():
        k_k2ydata = model_get_ydata2(objs, attr)
        for k, k2ydata in k_k2ydata.iteritems():
            if not isinstance(k, str):
                k = str(k)
            outbase = os.path.join(outdir, "%s_%s" % (attr, k))
            if attr in ['v2del', 'j2del', 'd2del']:
                rcmp.diff_plot(k2ydata, outbase, label=attr.upper(), xmin=-0.5,
                          xmax=16.5, ylabel='Frequency')
            else:
                rcmp.diff_plot(k2ydata, outbase, label=attr.upper(),
                               ylabel='Frequency')

def main():
    indir = sys.argv[1]
    outdir = sys.argv[2]
    
    model_makeplots(indir, outdir)

if __name__ == '__main__':
    main()


