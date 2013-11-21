#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a TCR repertoire sample
'''

import random
import copy


class Sample():
    '''Represents a sample
    '''
    def __init__(self, name, clones=None, group=None):
        self.name = name
        if clones is None:
            self.clones = []
        else:
            self.clones = clones
        self.group = group
        self.size = 0  # total sequence (read) count
        self.size = sum([clone.count for clone in self.clones])
        self.numclone = len(self.clones)

    def addclone(self, clone):
        if clone is not None:
            self.clones.append(clone)
            self.size = self.size + clone.count
            self.numclone = self.numclone + 1

    def addclones(self, clones):
        for clone in clones:
            self.addclone(clone)

    def setgroup(self, group):
        self.group = group

    def setclones(self, clones):
        self.clones = []
        self.size = 0
        self.numclone = 0
        if clones is not None:
            self.addclones(clones)

    def resetfreq(self):
        if len(self.clones) > 0:
            self.numclone = len(self.clones)
            self.size = sum([clone.count for clone in self.clones])
            for clone in self.clones:
                clone.freq = float(clone.count)/self.size

#======= RELEVANT FUNCTIONS =======
def filter_by_size(sample, minsize, maxsize=-1, freq=False, freqadjust=False):
    # Remove clones that have count (freq) < minsize 
    # or count (freq) > maxsize if maxsize is specified (!= -1)
    # if freqadjust=True, recompute clones' frequencies in the filtered
    # repertoire
    newsample = Sample(sample.name, group=sample.group)
    newclones = []
    for clone in sample.clones:
        newclone = copy.copy(clone)
        size = newclone.count
        if freq:
            size = newclone.freq
        if size >= minsize:
            if maxsize == -1 or size <= maxsize:
                newclones.append(newclone)
    newsample.setclones(newclones)
    if freqadjust:
        newsample.resetfreq()
    return newsample

#======= SAMPLING =============
'''
Different samples have different level of sequencing, and hence will create
bias. For example, large sample (lots of sequences got amplified and sequenced)
will inherently have more overlapping with other samples than a smaller sample.
The purpose of sampling is to make sure that every sample has the same number 
of starting sequences to avoid the bias
'''
def sampling(sample, size):
    if sample is None or sample.size == 0: 
        raise ValueError("Sample sampling: sample is None or has 0 clone.\n")
    if size <= 0 or size > sample.size:
        raise ValueError("Sample sampling: invalid size %d.\n" % size)
    
    newsample = Sample(sample.name, group=sample.group)
    indices = []
    for i, clone in enumerate(sample.clones):
        indices.extend([i] * clone.count)
    chosen_indices = random.sample(indices, size)

    index2count = {}
    for i in chosen_indices:
        if i not in index2count:
            index2count[i] = 1
        else:
            index2count[i] = index2count[i] + 1
    
    newclones = []
    for i, count in index2count.iteritems():
        newclone = copy.copy(sample.clones[i])
        newclone.count = count
        newclone.freq = float(count)/size
        newclones.append(newclone)
    newsample.setclones(newclones)

    return newsample

def sampling_uniq(sample, size):
    # randomly select "size" number of clones. 
    # Note: Use with caution -- each clone's count and freq stays the
    # same (i.e sum of freq may not add up to 1)
    if sample is None or sample.size == 0:
        raise ValueError("Sample sampling: sample is None or has 0 clone.\n")
    if size <= 0 or size > len(sample.clones):
        raise ValueError("Sample sampling: invalid size %d.\n" % size)

    newsample = Sample(sample.name, group=sample.group)
    newclones = random.sample(sample.clones, size)
    newsample.setclones(newclones)

    return newsample
    
        

