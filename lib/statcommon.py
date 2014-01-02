#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common stats functions
'''

import os
import sys

from scipy.stats import ttest_ind
from scipy.stats import ttest_rel
from scipy.stats import fisher_exact
import numpy as np

class SampleStat:
    def __init__(self):
        #sample info
        self.name = None
        self.group = None
        self.color = None
        self.marker = None
        self.size = 0
        self.numclone = 0
        
    def set_sample_info(self, sample):
        self.name = sample.name
        self.group = sample.group
        self.color = sample.color
        self.marker = sample.marker
        self.size = sample.size
        self.numclone = sample.numclone

    def set_sample_info2(self, name, group, color, marker, size=None,
                                                                numclone=None):
        self.name = name
        self.group = group
        self.color = color
        self.marker = marker
        if size is not None:
            self.size = size
        if numclone is not None:
            self.numclone = numclone
        
    def getitems(self):
        return self.__dict__.keys()
    
    def __getitem__(self, name):
        if name not in self.__dict__:
            return None
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

class PairStat:
    '''Stats of a pair of samples
    '''
    def __init__(self):
        #sample info
        self.name1 = None  # name of sample 1
        self.name2 = None  
        self.group1 = None  # group that sample 1 belongs to
        self.group2 = None
        self.color = None
        self.marker = None
        
    def set_sample_info(self, sample1, sample2):
        self.name1 = sample1.name
        self.group1 = sample1.group
        self.name2 = sample2.name
        self.group2 = sample2.group
        if self.group1 == self.group2:
            self.color = "#FE8E8F"
            self.marker = "o"
        else:
            self.color = "#A6D7FE"
            self.marker = "^"

    def set_sample_info2(self, name1, group1, name2, group2, color, marker):
        self.name1 = name1
        self.group1 = group1
        self.name2 = name2
        self.group2 = group2
        self.color = color
        self.marker = marker
        
    def getitems(self):
        return self.__dict__.keys()
    
    def __getitem__(self, name):
        if name not in self.__dict__:
            return None
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

#========== FUNCTIONS ==========
def get_clone2samples(samples, sizetype='count'):
    clone2sample2size = {}
    for sample in samples:
        name = sample.name
        for clone in sample.clones:
            ids = clone.get_vjseq_ids()
            assert len(ids) > 0
            if sizetype not in clone.getitems():
                raise KeyError("Clone does not have attribute %s" % sizetype)
            size = clone[sizetype]
            size /= len(ids)
            if size == 0 and sizetype == 'count':  #len(ids) > size
                size = 1
            for id in ids:
                if id not in clone2sample2size:
                    clone2sample2size[id] = {name: size}
                elif name not in clone2sample2size[id]:
                    clone2sample2size[id][name] = size
                else:
                    clone2sample2size[id][name] += size
    return clone2sample2size

def sizedict_to_vec(mydict, freq=False):
    # flatten val2size dictionary (key=a value, val=count of that val)
    # to a vector. Ex: {'a':2, 'b':1} --> ['a', 'a', 'b']
    vec = []
    for k, v in mydict.iteritems():
        if freq:
            v = int(v * 100)
        else:
            v = int(v)
        vec.extend([k] * v)
    return vec

def obj_dictattr_lookup(obj, args):
    # obj has a dictionary attribute
    assert len(args) == 2
    attr = args[0]
    key = args[1]
    if key in obj[attr]:
        return obj[attr][key]
    else:
        return None

def group_vector(names, name2obj, attr=None, func=None, func_args=None):
    # args is optional, additional arguments
    vec = []
    for n in names:
        obj = name2obj[n]
        if attr:
            vec.append(obj[attr])
        elif func:
            vec.append(func(obj, func_args))
        else:
            vec.append(obj)
    return vec

def pair_matrix(rnames, cnames, pair2obj, attr=None, sep="_", func=None, args=None):
    # get list of list (matrix) from pairwise comparisons
    # sep = separator/delimitor
    rows = []
    for n1 in rnames:  # row
        row = []
        for n2 in cnames:  # column
            pair = "%s%s%s" % (n1, sep, n2)
            pairrv = "%s%s%s" % (n2, sep, n1)
            obj = None
            if pair in pair2obj:
                obj = pair2obj[pair]
            elif pairrv in pair2obj:
                obj = pair2obj[pairrv]
            if obj:
                if attr:
                    row.append(obj[attr])
                elif func:
                    row.append(func(obj, args))
                else:
                    row.append(obj)
            else:
                row.append(0.0)
        rows.append(row)
    return rows

def pair_vec(names1, names2, pair2stat, attr):
    # return vector of all pairwise comparisition from group names1
    # and group names2
    vec = []
    visited = []
    for n1 in names1:
        for n2 in names2:
            if n1 == n2:
                continue
            pair = "%s_%s" % (n1, n2)
            pair_rv = "%s_%s" % (n2, n1)
            if pair in visited or pair_rv in visited:
                continue
            obj = None
            if pair in pair2stat:
                obj = pair2stat[pair]
            elif pair_rv in pair2stat[pair_rv]:
                obj = pair2stat[pair_rv]
            if obj:
                vec.append(obj[attr])
                visited.append(pair)
    return vec
    
def ttest_allpairs(group2names, name2obj, matched, attr=None, func=None,
                   func_args=None): 
    # perform ttests for all pairs of groups
    assert group2names and len(group2names) >= 2
    pair2stats = {}  # key=group1_group2; val=(t, p)
    group2mean = {}  # key=group; val=(mean, std)
    ttest = ttest_ind
    if matched:
        ttest = ttest_rel
    groups = group2names.keys()
    for i1 in xrange(0, len(groups) -1):
        g1 = groups[i1]
        vec1 = group_vector(group2names[g1], name2obj, attr, func, func_args)
        group2mean[g1] = (np.mean(vec1), np.std(vec1))
        for i2 in xrange(i1 + 1, len(groups)):
            g2 = groups[i2]
            vec2 = group_vector(group2names[g2], name2obj, attr, func,
                                                                     func_args)
            if i2 == len(groups) - 1:
                group2mean[g2] = (np.mean(vec2), np.std(vec2))
            pair = "%s_%s" % (g1, g2)
            tval, pval = ttest(vec1, vec2)
            pair2stats[pair] = (tval, pval)
    return pair2stats, group2mean

def ttest_pair(vec1, vec2, matched=False):
    ttest = ttest_ind
    if matched:
        assert len(vec1) == len(vec2)
        ttest = ttest_rel
    tval, pval = ttest(vec1, vec2)
    return tval, pval

def ttest_write(f, name, pair2tp, group2mean, pcutoff=1):
    for pair, (tval, pval) in pair2tp.iteritems():
        if pval <= pcutoff:
            f.write("%s\t%s\t%.2e\t%.2e" % (name, pair, tval, pval))
            groups = pair.split("_")
            (m1, std1) = group2mean[groups[0]]
            (m2, std2) = group2mean[groups[1]]
            f.write("\t%.2e +/- %.2e\t%.2e +/- %.2e\n" % (m1, std1, m2, std2))
        
def vec_mean(vec):
    return np.mean(vec)

def vec_std(vec):
    return np.std(vec)

def vec_mean_std(vec):
    return vec_mean(vec), vec_std(vec)

#def fisher_exact_test():
#    # cols = <Group1> <Group2>; rows = <Present> <Absent>

