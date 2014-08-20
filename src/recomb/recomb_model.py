#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''Compute recombination probabilities
E.g: for each sample repertoire, compute P(V), P(J), P(D, J),
P(delV|V), P(delJ|J), P(del5, del3|D), P(insVD), P(x_i|x_i-1),
P(insDJ), P(y_j|y_j+1)
Compute the median of a subset of samples etc...
'''

import os
import sys
import gzip
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from sonLib.bioio import system

import aimseqtk.lib.sample as libsample
#import aimseqtk.src.recomb.recomb_common as rcommon
from aimseqtk.src.recomb.recomb_common import RecombModel


####################### Functions ############################
def dict_convert_to_freq(mydict):
    total = sum(mydict.values())
    if total > 0:
        for k, v in mydict.iteritems():
            mydict[k] = float(v) / total

def dict_convert_to_freq_depth2(mydict):
    for k, dict2 in mydict.iteritems():
        dict_convert_to_freq(dict2)

def model_convert_to_freq(model):
    for attr in ['v', 'j', 'd', 'dj', 'ins_vd', 'ins_dj']:
        dict_convert_to_freq(model[attr])
    for attr in ['v2del', 'j2del', 'd2del', 'vd2next', 'dj2prev']:
        dict_convert_to_freq_depth2(model[attr]) 

def update_dict(dict1, dict2):
    # dictionary has 1 level: key 2 value
    for k, v in dict2.iteritems():
        if k not in dict1:
            dict1[k] = v
        else:
            dict1[k] += v

def update_dict_depth2(dict1, dict2):
    # dictionary has 2 depth levels: {k1: {k2: v}}
    for k1, k2v in dict2.iteritems():
        if k1 not in dict1:
            dict1[k1] = k2v
        else:
            update_dict(dict1[k1], k2v) 

def model_update(model1, model2):
    for attr in ['v', 'j', 'd', 'dj', 'ins_vd', 'ins_dj']:
        update_dict(model1[attr], model2[attr])
    for attr in ['v2del', 'j2del', 'd2del', 'vd2next', 'dj2prev']:
        update_dict_depth2(model1[attr], model2[attr])

def get_recomb_stats0(clones):
    # in case of unambiguous calls, just record the first call
    model = RecombModel()
    for clone in clones:
        # V
        model.update_attr('v', clone.vgenes[0], 1)
        model.update_attr2('v2del', clone.vgenes[0], clone.vdel, 1)
        # J
        model.update_attr('j', clone.jgenes[0], 1)
        model.update_attr2('j2del', clone.jgenes[0], clone.jdel, 1)
        # D
        model.update_attr('d', clone.dgenes[0], 1)
        model.update_attr2('d2del', clone.dgenes[0],
                           (clone.d5del, clone.d3del), 1)
        # DJ
        model.update_attr('dj', (clone.dgenes[0], clone.jgenes[0]), 1)
        # Insertion length
        ins_vd = clone.firstdpos - clone.lastvpos - 1
        model.update_attr('ins_vd', ins_vd, 1)
        ins_dj = clone.firstjpos - clone.lastdpos - 1
        model.update_attr('ins_dj', ins_dj, 1)
        # inserted nt given 5' nt
        vd_nts = clone.nuc[clone.lastvpos: clone.firstdpos]  # include the lastV
        model.update_vd2next(vd_nts, 1)
        dj_nts = clone.nuc[clone.lastdpos + 1: clone.firstjpos + 1]  # include lastJ
        model.update_dj2prev(dj_nts, 1)
    return model

def get_recomb_stats_splitweight0(clones):
    # in case of unambiguous calls, split the weight
    model = RecombModel()
    for clone in clones:
        # V
        vsize = 1.0 / len(clone.vgenes)
        for v in clone.vgenes:
            model.update_attr('v', v, vsize)
        model.update_attr2('v2del', clone.vgenes[0], clone.vdel, vsize)
        # J
        jsize = 1.0 / len(clone.jgenes)
        for j in clone.jgenes:
            model.update_attr('j', j, jsize)
        model.update_attr2('j2del', clone.jgenes[0], clone.jdel, jsize)
        # D
        dsize = 1.0 / len(clone.dgenes)
        for d in clone.dgenes:
            model.update_attr('d', d, dsize)
        model.update_attr2('d2del', clone.dgenes[0], (clone.d5del, clone.d3del), dsize)
        # DJ
        numdj = len(clone.jgenes) * len(clone.dgenes)
        djsize = 1.0 / numdj
        for d in clone.dgenes:
            for j in clone.jgenes:
                model.update_attr('dj', (d, j), djsize)
        # Insertion length
        ins_vd = clone.firstdpos - clone.lastvpos - 1
        model.update_attr('ins_vd', ins_vd, vsize)
        ins_dj = clone.firstjpos - clone.lastdpos - 1
        model.update_attr('ins_dj', ins_dj, jsize)
        # inserted nt given 5' nt
        vd_nts = clone.nuc[clone.lastvpos: clone.firstdpos]  # include the lastV
        model.update_vd2next(vd_nts, vsize)
        dj_nts = clone.nuc[clone.lastdpos + 1: clone.firstjpos + 1]  # include lastJ
        model.update_dj2prev(dj_nts, jsize)
    return model

def get_recomb_stats(clones):
    # in case of unambiguous calls, just record the first call
    model = RecombModel()
    for clone in clones:
        # V
        model.update_attr('v', clone.v, 1)
        model.update_attr2('v2del', clone.v, clone.vdel, 1)
        # J
        model.update_attr('j', clone.j, 1)
        model.update_attr2('j2del', clone.j, clone.jdel, 1)
        # D
        model.update_attr('d', clone.d, 1)
        model.update_attr2('d2del', clone.d, (clone.d5del, clone.d3del), 1)
        # DJ
        model.update_attr('dj', (clone.d, clone.j), 1)
        # Insertion length
        ins_vd = len(clone.vdins) - 1
        model.update_attr('ins_vd', ins_vd, 1)
        ins_dj = len(clone.djins) - 1
        model.update_attr('ins_dj', ins_dj, 1)
        # inserted nt given 5' nt
        model.update_vd2next(clone.vdins, 1)  # include the lastV
        model.update_dj2prev(clone.djins, 1)  # include the fistJ
    return model

def get_recomb_stats_splitweight(clones):
    # in case of unambiguous calls, split the weight
    model = RecombModel()
    for clone in clones:
        # V
        model.update_attr('v', clone.v, 1.0)
        model.update_attr2('v2del', clone.v, clone.vdel, 1.0)
        # J
        model.update_attr('j', clone.j, 1.0)
        model.update_attr2('j2del', clone.j, clone.jdel, 1.0)
        # D
        model.update_attr('d', clone.d, 1.0)
        model.update_attr2('d2del', clone.d, (clone.d5del, clone.d3del), 1.0)
        # DJ
        model.update_attr('dj', (clone.d, clone.j), 1.0)
        # Insertion length
        ins_vd = len(clone.vdins) - 1
        model.update_attr('ins_vd', ins_vd, 1.0)
        ins_dj = len(clone.djins) - 1
        model.update_attr('ins_dj', ins_dj, 1.0)
        # inserted nt given 5' nt
        model.update_vd2next(clone.vdins, 1.0)
        model.update_dj2prev(clone.djins, 1.0)
    return model

def write_attr(mydict, outfile):
    f = open(outfile, 'w')
    keys = sorted(mydict.keys())
    f.write("#Name\tFreq\n")
    for k in keys:
        f.write("%s\t%f\n" % (k, mydict[k]))
    f.close()

def write_attr2(mydict, outfile):
    f = open(outfile, 'w')
    cols = sorted(mydict.keys())
    f.write("#\t%s\n" % "\t".join(cols))
    rows = []
    for col in cols:
        for r in mydict[col].keys():
            if r not in rows:
                rows.append(r)
    for row in sorted(rows):
        f.write("%s" % str(row))
        for col in cols:
            if row in mydict[col]:
                f.write("\t%f" % mydict[col][row])
            else:
                f.write("\t0.0")
        f.write("\n")
    f.close()
    
def model_write(model, outdir):
    # write model probabilites to text files
    for attr in ['v', 'j', 'd', 'dj', 'ins_vd', 'ins_dj']:
        outfile = os.path.join(outdir, "%s.txt" % attr)
        write_attr(model[attr], outfile)
    for attr in ['v2del', 'j2del', 'd2del', 'vd2next', 'dj2prev']:
        outfile = os.path.join(outdir, "%s.txt" % attr)
        write_attr2(model[attr], outfile)
      
class ClonesRecombModel(Target):
    '''Get the recomb_model related counts for a subset of clones
    Return a "RecomModel" obj picked to outfile
    '''
    def __init__(self, infile, outfile):
        Target.__init__(self)
        self.infile = infile
        self.outfile = outfile

    def run(self):
        clones = pickle.load(gzip.open(self.infile, 'rb'))
        recomb_model = get_recomb_stats(clones)
        pickle.dump(recomb_model, gzip.open(self.outfile, 'wb'))

class SampleRecombModelAgg(Target):
    '''Combine stats of each batched computed indepently to 1 model
    Convert count to frequencies (probabilities)
    '''
    def __init__(self, indir, outdir):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        model = RecombModel()
        for file in os.listdir(self.indir):
            filepath = os.path.join(self.indir, file)
            submodel = pickle.load(gzip.open(filepath, 'rb'))
            model_update(model, submodel)
        # convert to frequecies
        model_convert_to_freq(model)
        outfile = os.path.join(self.outdir, "model.pickle")
        pickle.dump(model, gzip.open(outfile, "wb"))
        model_write(model, self.outdir)

class SampleRecombModel(Target):
    '''Get the recombination model for the input sample
    '''
    def __init__(self, indir, outdir):
        Target.__init__(self)
        self.indir = indir  # indir contains pickle files (batches) of clones 
        self.outdir = outdir  # directory to put outfile file there

    def run(self):
        name = os.path.basename(self.indir.rstrip("/"))
        global_dir = self.getGlobalTempDir()
        tempdir = os.path.join(global_dir, "recom_model", name)
        system("mkdir -p %s" % tempdir)

        for file in os.listdir(self.indir):  # each batch
            if file == os.path.basename(self.indir):
                continue
            infile = os.path.join(self.indir, file)
            outfile = os.path.join(tempdir, file)
            self.addChildTarget(ClonesRecombModel(infile, outfile))
        self.setFollowOnTarget(SampleRecombModelAgg(tempdir, self.outdir))








