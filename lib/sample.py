#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a TCR repertoire sample
'''

import random
import copy
import time
import gzip
import cPickle as pickle

from jobTree.scriptTree.target import Target
from sonLib.bioio import logger


class Sample():
    '''Represents a sample
    '''
    def __init__(self, name, clones=None, group=None, color=(0, 0, 0), 
                                                           marker='.'):
        self.name = name
        if clones is None:
            self.clones = []
        else:
            self.clones = clones
        self.group = group
        self.size = 0  # total sequence (read) count
        self.size = sum([clone.count for clone in self.clones])
        self.numclone = len(self.clones)
        self.color = color
        self.marker = marker

    def __lt__(self, other):
        return self.size < other.size
    
    def __le__(self, other):
        return self.size <= other.size

    def __eq__(self, other):
        return self.size == other.size

    def __ne__(self, other):
        return self.size != other.size

    def __gt__(self, other):
        return self.size > other.size

    def __ge__(self, other):
        return self.size >= other.size

    def getitems(self):
        return self.__dict__.keys()

    def __getitem__(self, name):
        if name not in self.__dict__:
            return None
        return self.__dict__[name]
    
    def __setitem__(self, name, val):
        self.__dict__[name] = val

    def addclone(self, clone):
        if clone is not None:
            self.clones.append(clone)
            self.size = self.size + clone.count
            self.numclone = self.numclone + 1

    def addclones(self, clones):
        count = sum([clone.count for clone in clones])
        self.clones.extend(clones)
        self.size += count
        self.numclone += len(clones)

    def setgroup(self, group):
        self.group = group

    def setcolor(self, color):
        self.color = color

    def setmarker(self, marker):
        self.marker = marker

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

######## RELEVANT FUNCTIONS ########
#======== Write a sample to output file =======
def write_sample(file, sample):
    if not sample.clones or sample.numclone <= 0 or sample.size <= 0:
        return
    assert len(sample.clones) > 0

    f = open(file, 'w')
    firstclone = sample.clones[0]
    columns = firstclone.get_sorted_items()
    f.write("%s\n" % "\t".join(columns))
    for clone in sample.clones:
        f.write("%s\n" % clone.getstr())
    f.close()

def write_samples(outdir, samples):
    for sample in samples:
        outfile = os.path.join(outdir, "%s.tsv" % sample.name)
        write_sample(outfile, sample)

#======= Filtering =============
def filter_by_size(sample, mincount=-1, maxcount=-1, minfreq=-1, maxfreq=-1,
                   freqadjust=False):
    # Remove clones that have count (freq) < minsize 
    # or count (freq) > maxsize if maxsize is specified (!= -1)
    # if freqadjust=True, recompute clones' frequencies in the filtered
    # repertoire
    newsample = Sample(sample.name, group=sample.group, color=sample.color,
                                                        marker=sample.marker)
    newclones = []
    for clone in sample.clones:
        newclone = copy.deepcopy(clone)
        if newclone.count >= mincount and newclone.freq >= minfreq:
            if ((maxcount == -1 or newclone.count <= maxcount) and
                (maxfreq == -1 or newclone.freq <= maxfreq)):
                newclones.append(newclone)
    newsample.setclones(newclones)
    if freqadjust:
        newsample.resetfreq()
    return newsample

def filter_by_status(sample, productive=True, resetfreq=True):
    newsample = Sample(sample.name, group=sample.group, color=sample.color,
                                                        marker=sample.marker)
    newclones = []
    for clone in sample.clones:
        if productive == clone.productive:
            newclones.append(copy.deepcopy(clone))
    newsample.setclones(newclones)
    if resetfreq:
        newsample.resetfreq()
    return newsample

class WriteSample(Target):
    def __init__(self, psample, npsample, samout, outdir):
        self.psample = psample
        self.npsample = npsample
        self.samout = samout
        self.outdir = outdir

    def run(self):
        pdir = os.path.join(self.outdir, "samples", "productive")
        npdir = os.path.join(self.outdir, "samples", "non_productive")
        name = self.psample.name
        if 'pickle' in self.samout:
            pfile = os.path.join(pdir, "pickle", "%s.pickle" % name)
            pickle.dump(self.psample, gzip.open(pfile, "wb"))
            npfile = os.path.join(npdir, "pickle", "%s.pickle" % name)
            pickle.dump(self.npsample, gzip.open(pfile, "wb"))
        if 'txt' in self.samout:
            pfile = os.path.join(pdir, "txt", "%s.txt" % name)
            write_sample(pfile, self.psample) 
            npfile = os.path.join(npdir, "txt", "%s.txt" % name)
            write_sample(npfile, self.npsample) 
                
class FilterBySize(Target):
    def __init__(self, outfile, sample, mincount=-1, maxcount=-1,
                 minfreq=-1, maxfreq=-1, freqadjust=False):
        Target.__init__(self)
        self.outfile = outfile
        self.sample = sample
        self.mincount = mincount
        self.maxcount = maxcount
        self.minfreq = minfreq
        self.maxfreq = maxfreq
        self.freqadjust = freqadjust

    def run(self):
        starttime = time.time()
        filterbysize = False
        if (self.mincount > 1 or self.maxcount > 0 or self.minfreq > 0 or
            self.maxfreq > 0):
            filterbysize = True
        sample = self.sample
        if filterbysize:
            sample = filter_by_size(self.sample, self.mincount, self.maxcount, 
                                  self.minfreq, self.maxfreq, self.freqadjust)
        pickle.dump(sample, gzip.open(self.outfile, "wb")) 
        mytime = time.time() - starttime
        logger.debug("Filter_by_size for sample %s in %.4f seconds" % 
                     (self.sample.name, mytime))
    
class FilterByStatus(Target):
    def __init__(self, pdir, npdir, sample, resetfreq=True, options=None):
        Target.__init__(self)
        self.pdir = pdir
        self.npdir = npdir
        self.sample = sample
        self.resetfreq = resetfreq
        self.options = options

    def run(self):
        starttime = time.time()
        productive_sample = filter_by_status(self.sample, True, self.resetfreq)
        nonproductive_sample = filter_by_status(self.sample, False, 
                                                self.resetfreq)
        pfile = os.path.join(self.pdir, "%s.pickle" % self.sample.name)
        pickle.dump(productive_sample, gzip.open(pfile, "wb"))
        npfile = os.path.join(self.npdir, "%s.pickle" % self.sample.name)
        pickle.dump(nonproductive_sample, gzip.open(npfile, "wb"))
        
        # Write to outdir if requested
        opts = self.options
        if opts and opts.samout:
            self.addChildTarget(WriteSample(productive_sample,
                            nonproductive_sample, opts.samout, opts.outdir))
        mytime = time.time() - starttime
        logger.debug("Filter_by_status for sample %s in %.4f seconds" %
                     (self.sample.name, mytime))

#======= set group, color, marker etc =====
def set_sample_group(sample, name2group):
    if sample.name in name2group:
        sample.setgroup(name2group[sample.name])

def set_sample_color(sample, name2color):
    if sample.name in name2color:
        sample.setcolor(name2color[sample.name])

def set_sample_marker(sample, name2marker):
    if sample.name in name2marker:
        sample.setmarker(name2marker[sample.name])

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
    
    newsample = Sample(sample.name, group=sample.group, color=sample.color,
                                                        marker=sample.marker)
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
        newclone = copy.deepcopy(sample.clones[i])
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

    newsample = Sample(sample.name, group=sample.group, color=sample.color,
                                                        marker=sample.marker)
    newclones = random.sample(sample.clones, size)
    newsample.setclones(newclones)

    return newsample

class Sampling(Target):
    def __init__(self, sample, size, outdir):
        Target.__init__(self)
        self.sample = sample
        self.size = size
        self.outdir = outdir

    def run(self):
        subsample = sampling(self.sample, self.size)
        outfile = os.path.join(self.outdir, "%s.pickle" % sample.name)
        pickle.dump(subsample, gzip.open(outfile, "wb"))
       

