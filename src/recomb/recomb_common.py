
import sys
from math import log10
import numpy as np
import aimseqtk.lib.common as lcommon


class RecombModel:
    '''Represents recombination/ generative model of a repertoire
    '''
    def __init__(self):
        self.v = {}  # key=specific v, val=frequencies/count
        self.j = {}
        self.d = {}
        self.dj = {}
        self.v2del = {}  # delV_length | V
        self.j2del = {}  # delJ_length | J
        self.d2del = {}  # delD5_delD3 | D
        self.ins_vd = {}  # insertion length btw V and D
        self.vd2next = {}  # {5'_VD_nucleotide: {VD_nucleotide: freq}}
        self.ins_dj = {}  # insertion length btw D and J
        self.dj2prev = {}
    
    def get_attrs_depth1(self):
        return ['v', 'd', 'j', 'dj', 'ins_vd', 'ins_dj']

    def get_attrs_depth2(self):
        return ['v2del', 'j2del', 'd2del', 'vd2next', 'dj2prev']
    
    def getitems(self):
        return self.__dict__.keys()

    def __getitem__(self, name):
        if name not in self.__dict__:
            return None
        return self.__dict__[name]
    
    def __setitem__(self, name, val):
        self.__dict__[name] = val

    def update_vd2next(self, nts, size):
        if not nts or len(nts) < 2:
            return
        for i in xrange(1, len(nts)):
            prev_nt = nts[i - 1]
            nt = nts[i]
            if prev_nt not in self.vd2next:
                self.vd2next[prev_nt] = {nt: size}
            elif nt not in self.vd2next[prev_nt]:
                self.vd2next[prev_nt][nt] = size
            else:
                self.vd2next[prev_nt][nt] += size

    def update_dj2prev(self, nts, size):
        if not nts or len(nts) < 2:
            return
        for i in xrange(len(nts) - 1, 0, -1):
            nt = nts[i]
            prev_nt = nts[i - 1]
            if nt not in self.dj2prev:
                self.dj2prev[nt] = {prev_nt: size}
            elif prev_nt not in self.dj2prev[nt]:
                self.dj2prev[nt][prev_nt] = size
            else:
                self.dj2prev[nt][prev_nt] += size

    def update_attr(self, attr, key, size):
        assert attr in ['v', 'd', 'j', 'dj', 'ins_vd', 'ins_dj']
        if key not in self[attr]:
            self[attr][key] = size
        else:
            self[attr][key] += size

    def update_attr2(self, attr, key1, key2, size):
        assert attr in ['v2del', 'j2del', 'd2del']
        #if isinstance(key2, int):
        if key1 not in self[attr]:
            self[attr][key1] = {key2: size}
        elif key2 not in self[attr][key1]:
            self[attr][key1][key2] = size
        else:
            self[attr][key1][key2] += size

def union_models(models):
    # for each attribute, find the union list of keys and add the
    # missing keys of to each sample (with val = 0)
    if not models:
        return
    for attr in models[0].get_attrs_depth1():
        key_lists = [m[attr].keys() for m in models]
        union_keys = lcommon.union_lists(key_lists)
        for k in union_keys:
            for m in models:
                if k not in m[attr]:
                    m[attr][k] = 0.0
    for attr in models[0].get_attrs_depth2():
        key_lists = [m[attr].keys() for m in models]
        union_keys = lcommon.union_lists(key_lists)
        for k in union_keys:
            for m in models:
                if k not in m[attr]:
                    m[attr][k] = {}
        for k in union_keys:
            key_lists2 = [m[attr][k].keys() for m in models]
            union_keys2 = lcommon.union_lists(key_lists2)
            for k2 in union_keys2:
                for m in models:
                    if k2 not in m[attr][k]:
                        m[attr][k][k2] = 0.0

def model_get_median(models):
    # note: has to do "union_models" before this
    if not models:
        return None
    elif len(models) == 1:
        return models[0]

    med_model = RecombModel()
    for attr in med_model.get_attrs_depth1():
        for k in models[0][attr].keys():
            vec = [m[attr][k] for m in models]
            med_model[attr][k] = np.median(vec)

    for attr in med_model.get_attrs_depth2():
        for k, mydict in models[0][attr].iteritems():
            med_model[attr][k] = {}
            for k2 in mydict.keys():
                vec = [m[attr][k][k2] for m in models]
                med_model[attr][k][k2] = np.median(vec)

    return med_model

def get_median_model(model_dir):
    objs = lcommon.load_pickledir(model_dir)
    union_models(objs)
    return model_get_median(objs)

def get_log10(v):
    if v == 0:
        return float('-inf')
    else:
        return log10(v)

def ntclone_likelihood(ntclone, model):
    try:
        #if ntclone.v not in model.v:
        #    sys.stderr.write("V %s is not in model\n" % ntclone.v)
        #    sys.exit()
        v = model.v[ntclone.v]
        #if (ntclone.d, ntclone.j) not in model.dj:
        #    sys.stderr.write("(%s, %s) is not in model.dj\n" % (ntclone.d, ntclone.j))
        #    print model.dj
        #    sys.exit()
        dj = model.dj[(ntclone.d, ntclone.j)]
        #if ntclone.v not in model.v2del:
        #    sys.stderr.write("V %s is not in model\n" % (ntclone.v))
        #    sys.exit()
        #elif ntclone.vdel not in model.v2del[ntclone.v]:
        #    sys.stderr.write("Vdel %d not in v2del[%s]\n" % (ntclone.vdel, ntclone.v))
        #    sys.exit()
        vdel = model.v2del[ntclone.v][ntclone.vdel]
        #if ntclone.j not in model.j2del:
        #    sys.stderr.write("V %s is not in model\n" % (ntclone.j))
        #    sys.exit()
        #elif ntclone.jdel not in model.j2del[ntclone.j]:
        #    sys.stderr.write("Vdel %d not in j2del[%s]\n" % (ntclone.jdel, ntclone.j))
        #    sys.exit()
        jdel = model.j2del[ntclone.j][ntclone.jdel]
        ddel = 1.0
        if ntclone.d in model.d2del:
            #if (ntclone.d5del, ntclone.d3del) not in model.d2del[ntclone.d]:
            #    sys.stderr.write("(%d, %d) is not in d2del[%s]\n" % (ntclone.d5del, ntclone.d3del, ntclone.d))
            #    sys.exit()
            ddel = model.d2del[ntclone.d][(ntclone.d5del, ntclone.d3del)]
        
        vdins = ntclone.vdins
        if len(vdins) - 1 not in model.ins_vd:
            #sys.stderr.write("vdins %s has length %d that is not in model.ins_vd\n" % (vdins, len(vdins) - 1))
            #sys.exit()
            ins_vd = 0.0
        else:
            ins_vd = model.ins_vd[len(vdins) - 1]
            for i in xrange(1, len(vdins)):
                #assert vdins[i-1] in model.vd2next
                #assert vdins[i] in model.vd2next[vdins[i-1]]
                ins_vd *= model.vd2next[vdins[i - 1]][vdins[i]]
            
        djins = ntclone.djins
        if len(djins) - 1 not in model.ins_dj:
            #sys.stderr.write("djins %s has length %d that is not in model.ins_dj\n" % (djins, len(djins) - 1))
            #sys.exit()
            ins_dj = 0.0
        else:
            ins_dj = model.ins_dj[len(djins) - 1]
            for i in xrange(len(djins) - 1):
                #assert djins[i-1] in model.dj2prev
                #assert djins[i] in model.dj2prev[djins[i-1]]
                ins_dj *= model.dj2prev[djins[i + 1]][djins[i]]
        llhood = (get_log10(v) + get_log10(dj) + get_log10(vdel) + get_log10(jdel)
                  + get_log10(ddel) + get_log10(ins_vd) + get_log10(ins_dj))
        return llhood
    except:
        return float('-inf')

def cmp_clone(clone1, clone2):
    return (clone1.v == clone2.v and clone1.j == clone2.j and
            clone1.d == clone2.d and clone1.nuc == clone2.nuc and
            clone1.vdel == clone2.vdel and clone1.jdel == clone2.jdel and
            clone1.d5del == clone2.d5del and clone1.vdins == clone2.vdins and
            clone1.djins == clone2.djins)

def visited_event(events, clone):
    # return True if clone is already in "events"
    for e in events:
        if cmp_clone(e, clone):
            return True
    return False

