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
import numpy as np

def group_vector(names, name2obj, attr):
    vec = []
    for n in names:
        obj = name2obj[n]
        vec.append(obj[attr])
    return vec

def ttest_allpairs(group2names, name2obj, matched, attr): 
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
        vec1 = group_vector(group2names[g1], name2obj, attr)
        group2mean[g1] = (np.mean(vec1), np.std(vec1))
        for i2 in xrange(i1 + 1, len(groups)):
            g2 = groups[i2]
            vec2 = group_vector(group2names[g2], name2obj, attr)
            if i2 == len(groups) - 1:
                group2mean[g2] = (np.mean(vec2), np.std(vec2))
            pair = "%s_%s" % (g1, g2)
            tval, pval = ttest(vec1, vec2)
            pair2stats[pair] = (tval, pval)
    return pair2stats, group2mean

