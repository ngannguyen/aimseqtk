#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a TCR repertoire sample
'''

import os
import random
import copy
import time
import gzip
import cPickle as pickle

from jobTree.scriptTree.target import Target
from sonLib.bioio import system
from sonLib.bioio import logger

import aimseqtk.lib.clone as libclone
import aimseqtk.lib.common as libcommon


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

class EmptySampleError(Exception):
    pass

######## RELEVANT FUNCTIONS ########
#======== Write a sample to output file =======
def write_clones(file, clones, append=False):
    if not clones:
        return
    if append:
        f = open(file, 'a')
    else:
        f = open(file, 'w')
    firstclone = clones[0]
    columns = firstclone.get_sorted_items()
    f.write("%s\n" % "\t".join(columns))
    for clone in clones:
        f.write("%s\n" % clone.getstr())
    f.close()

#def write_samples(outdir, samples):
#    for sample in samples:
#        outfile = os.path.join(outdir, "%s.tsv" % sample.name)
#        write_sample(outfile, sample)

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

#======= Filtering =============
def filter_by_size(clones, mincount=-1, maxcount=-1, minfreq=-1, maxfreq=-1):
    # Remove clones that have count (freq) < minsize 
    # or count (freq) > maxsize if maxsize is specified (!= -1)
    
    newclones = []
    for clone in clones:
        newclone = copy.deepcopy(clone)
        if newclone.count >= mincount and newclone.freq >= minfreq:
            if ((maxcount == -1 or newclone.count <= maxcount) and
                (maxfreq == -1 or newclone.freq <= maxfreq)):
                newclones.append(newclone)
    return newclones

def filter_by_status(clones, productive=True):
    newclones = []
    for clone in clones:
        if productive == clone.productive:
            newclones.append(copy.deepcopy(clone))
    return newclones

#======= Split clones by VJ ===============
def split_clones_by_vj(clones, sample_name=None):
    v2j2clones = {}
    for clone in clones:
        numcombi = len(clone.vgenes) * len(clone.jgenes)
        if clone.dgenes:
            numcombi *= len(clone.dgenes)
        nuc = clone.cdr3nuc
        if nuc is None:
            nuc = clone.nuc
        count = clone.count/numcombi
        if count == 0:
            count = 1
        normcount = clone.normcount/numcombi
        freq = clone.freq/numcombi

        for v in clone.vgenes:
            for j in clone.jgenes:
                cdr3clones = []
                if not clone.dgenes:
                    cdr3clone = libclone.Cdr3Clone(count, nuc, v, j, '',
                                                   clone.cdr3aa, sample_name,
                                                   normcount, freq)
                    cdr3clones.append(cdr3clone)
                else:
                    for d in clone.dgenes:
                        cdr3clone = libclone.Cdr3Clone(count, nuc, v, j, d,
                                                    clone.cdr3aa, sample_name,
                                                    normcount, freq)
                        cdr3clones.append(cdr3clone)
                if v not in v2j2clones:
                    v2j2clones[v] = {j: cdr3clones}
                elif j not in v2j2clones[v]:
                    v2j2clones[v][j] = cdr3clones
                else:
                    v2j2clones[v][j].extend(cdr3clones)
    return v2j2clones
            
def sample_get_size(indir):
    name = os.path.basename(indir.rstrip('/'))
    samfile = os.path.join(indir, "%s" % name)
    sample = pickle.load(gzip.open(samfile, 'rb'))
    assert sample.name == name
    numclone = 0
    size = 0
    # get total size and total number of clones
    for vjfile in os.listdir(indir):
        if vjfile == name:
            continue
        clones = pickle.load(gzip.open(os.path.join(indir, vjfile), "rb"))
        numclone += len(clones)
        size += sum([clone.count for clone in clones])
    if size == 0:
        raise EmptySampleError("Sample %s at %s is empty" % (name, indir))
    sample.size = size
    sample.numclone = numclone
    pickle.dump(sample, gzip.open(samfile, "wb"))
    return numclone, size

def get_vj2sams(indir, sams=None):
    vj2sams = {}
    if sams is None:
        sams = os.listdir(indir)
    for sam in sams:
        samdir = os.path.join(indir, sam)
        for vj in os.listdir(samdir):
            if vj == sam:
                continue
            if vj not in vj2sams:
                vj2sams[vj] = [sam]
            else:
                vj2sams[vj].append(sam)
    return vj2sams

def sample_all_clones(samdir):
    name = os.path.basename(samdir.strip('/'))
    allclones = []
    for vj in os.listdir(samdir):
        if vj == name:
            continue
        vjfile = os.path.join(samdir, vj)
        clones = pickle.load(gzip.open(vjfile, 'rb'))
        allclones.extend(clones)
    return allclones

def reset_freqs_vj(infile, size):
    clones = pickle.load(gzip.open(infile, "rb"))
    for c in clones:
        c.freq = float(c.count)/size
    pickle.dump(clones, gzip.open(infile, "wb"))

#========== OBJs ====
class WriteSample(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        if os.path.exists(self.outfile):
            system("rm -f" % self.outfile)
        for batch in os.listdir(self.indir):
            batchfile = os.path.join(self.indir, batch)
            clones = pickle.load(gzip.open(batchfile, "rb"))
            write_clones(self.outfile, clones, True)

class WriteSamples(Target):
    def __init__(self, indir, outdir, samout):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.samout = samout

    def run(self):
        if 'pickle' in self.samout:
            pickledir = os.path.join(self.outdir, "pickle")
            system("mkdir -p %s" % pickledir)
            system("cp -r %s %s" % (self.indir, pickledir))
        if 'txt' in self.samout:
            txtdir = os.path.join(self.outdir, "txt")
            system("mkdir -p %s" % txtdir)
            for sam in os.listdir(self.indir):
                samdir = os.path.join(self.indir, sam)
                outfile = os.path.join(txtdir, sam)
                self.addChildTarget(WriteSample(samdir, outfile))

class FilterSample(Target):
    '''Filter a sample by clone size and by productive status
    '''
    def __init__(self, outdir, name, samplefile, opts):
        Target.__init__(self)
        self.outdir = outdir
        self.name = name
        self.samplefile = samplefile
        self.opts = opts

    def run(self):
        # filter by size
        starttime = time.time()
        opts = self.opts
        clones = pickle.load(gzip.open(self.samplefile, 'rb'))
        if (opts.mincount > 1 or opts.maxcount > 0 or opts.minfreq > 0 or
            opts.maxfreq > 0):
            clones = filter_by_size(clones, opts.mincount, opts.maxcount,
                                    opts.minfreq, opts.maxfreq)
        msg = ("Filter_by_size for file %s done in %.4f s" %
                                 (self.samplefile, time.time() - starttime))
        logger.info(msg)
        starttime = time.time()

        # filter by status
        pclones = filter_by_status(clones, True)
        npclones = filter_by_status(clones, False)
        
        filename = os.path.basename(self.samplefile)
        if pclones:
            pdir = os.path.join(self.outdir, "productive", self.name)
            system("mkdir -p %s" % pdir)
            pfile = os.path.join(pdir, filename)
            pickle.dump(pclones, gzip.open(pfile, "wb"))
        if npclones:    
            npdir = os.path.join(self.outdir, "non_productive", self.name)
            system("mkdir -p %s" % npdir)
            npfile = os.path.join(npdir, filename)
            pickle.dump(npclones, gzip.open(npfile, "wb"))
        msg = ("Filter_by_status for file %s done in %.4f s" %
                                 (self.samplefile, time.time() - starttime))
        logger.info(msg)
        self.setFollowOnTarget(libcommon.CleanupFile(self.samplefile))

class SplitClonesByV(Target):
    '''
    '''
    def __init__(self, infile, outdir, sample_name=None):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.sample_name = sample_name

    def run(self):
        clones = pickle.load(gzip.open(self.infile, "rb"))
        v2j2clones = split_clones_by_vj(clones, self.sample_name)
        for v, j2clones in v2j2clones.iteritems():
            vclones = []
            for j, vjclones in j2clones.iteritems():
                vclones.extend(vjclones)
            vfile = os.path.join(self.outdir, v)
            pickle.dump(vclones, gzip.open(vfile, "wb"))

class SplitClonesByVJ(Target):
    '''
    '''
    def __init__(self, infile, outdir, sample_name=None):
        Target.__init__(self)
        self.infile = infile
        self.outdir = outdir
        self.sample_name = sample_name

    def run(self):
        clones = pickle.load(gzip.open(self.infile, "rb"))
        v2j2clones = split_clones_by_vj(clones, self.sample_name)
        for v, j2clones in v2j2clones.iteritems():
            for j, vjclones in j2clones.iteritems():
                vjfile = os.path.join(self.outdir, "%s_%s" % (v, j))
                pickle.dump(vjclones, gzip.open(vjfile, "wb"))

class SampleResetFreqsVJ(Target):
    def __init__(self, infile, size):
        Target.__init__(self)
        self.infile = infile
        self.size = size

    def run(self):
        #stime = time.time()
        reset_freqs_vj(self.infile, self.size)
        #self.logToMaster("SampleResetFreqsVJ: done in %.4f s\n" %
        #                 (time.time() - stime))

class SampleResetFreqs(Target):
    '''
    '''
    def __init__(self, indir):
        Target.__init__(self)
        self.indir = indir

    def run(self):
        self.logToMaster("SampleResetFreqs\n")
        name = os.path.basename(self.indir.rstrip('/'))
        numclone, size = sample_get_size(self.indir)
        for vjfile in os.listdir(self.indir):
            if vjfile == name:
                continue
            vjpath = os.path.join(self.indir, vjfile)
            reset_freqs_vj(vjpath, size)
            #self.addChildTarget(SampleResetFreqsVJ(vjpath, size))

class VjSampleAgg(Target):
    def __init__(self, vj, batches, indir, outdir):
        Target.__init__(self)
        self.vj = vj
        self.batches = batches
        self.indir = indir
        self.outdir = outdir

    def run(self):
        #stime = time.time()
        vjfile = os.path.join(self.outdir, self.vj)
        clones = []
        for batch in self.batches:
            file = os.path.join(self.indir, batch, self.vj)
            currclones = pickle.load(gzip.open(file, "rb"))
            clones.extend(currclones)
        pickle.dump(clones, gzip.open(vjfile, "wb"))
        #self.logToMaster("VjSampleAgg: done in %.4f s\n" %
        #                 (time.time() - stime))
    
class MakeDbSampleAgg(Target):
    '''
    Each sample is stored in a directory, with:
    sample_name.pickle: contain Sample obj, without the clones
    V1J1.pickle: list of Cdr3Clones that have V1 and J1
    V1J2.pickle
    ...
    VnJm.pickle
    Also recompute the frequencies
    '''
    def __init__(self, indir, outdir, opts, group=None, color=(0, 0, 0),
                                                        marker='.'):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.opts = opts
        self.group = group
        self.color = color
        self.marker = marker

    def run(self):
        # create a Sample obj
        name = os.path.basename(self.indir.rstrip('/'))
        sample = Sample(name, group=self.group, color=self.color,
                                                marker=self.marker)
        samfile = os.path.join(self.outdir, name)
        pickle.dump(sample, gzip.open(samfile, 'wb'))
        
        # aggregate the batches
        vj2batches = get_vj2sams(self.indir)
        for vj, batches in vj2batches.iteritems():
            self.addChildTarget(VjSampleAgg(vj, batches, self.indir,
                                            self.outdir))
        
        self.setFollowOnTarget(SampleResetFreqs(self.outdir))
        
class MakeDbSample(Target):
    '''Set up smaller jobs to split clones by VJ
    Then set up follow_on job to merge them
    '''
    def __init__(self, indir, outdir, opts, group, color, marker):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.opts = opts
        self.group = group
        self.color = color
        self.marker = marker

    def run(self):
        name = os.path.basename(self.indir.rstrip("/"))
        
        global_dir = self.getGlobalTempDir()
        tempdir = os.path.join(global_dir, "split_by_vj", name)
        system("mkdir -p %s" % tempdir)
        
        for file in os.listdir(self.indir):  # each batch, split clone by VJ
            infile = os.path.join(self.indir, file)
            child_outdir = os.path.join(tempdir, os.path.splitext(file)[0])
            system("mkdir -p %s" % child_outdir)
            self.addChildTarget(SplitClonesByV(infile, child_outdir, name))
            #self.addChildTarget(SplitClonesByVJ(infile, child_outdir, name))
        self.setFollowOnTarget(MakeDbSampleAgg(tempdir, self.outdir, self.opts,
                                          self.group, self.color, self.marker)) 


#======= SAMPLING =============
'''
Different samples have different level of sequencing, and hence will create
bias. For example, large sample (lots of sequences got amplified and sequenced)
will inherently have more overlapping with other samples than a smaller sample.
The purpose of sampling is to make sure that every sample has the same number 
of starting sequences to avoid the bias
'''
class SamplingError(Exception):
    pass

def sampling(sample, sampledir, outdir, args=None):
    if not args:
        return sampledir
        #raise ValueError("Sample sampling: sample is None or has 0 clone.\n")
    size = args[0]

    if sample is None or sample.size == 0: 
        raise ValueError("Sample sampling: sample is None or has 0 clone.\n")
    if size <= 0 or size > sample.size:
        raise ValueError(("Sample sampling: invalid size %d.\n" % size +
                          "Sample %s has %d sequences.\n" % (sample.name, sample.size)))
    
    newsample = Sample(sample.name, group=sample.group, color=sample.color,
                                                        marker=sample.marker)
    indices = []  # represent all clones
    for vj in os.listdir(sampledir):
        if vj == sample.name:
            continue
        vjfile = os.path.join(sampledir, vj)
        clones = pickle.load(gzip.open(vjfile, 'rb'))
        for i, clone in enumerate(clones):
            indices.extend([(vj, i)] * clone.count)
    chosen_indices = random.sample(indices, size)

    vj2index2count = {}
    for (vj, i) in chosen_indices:
        if vj not in vj2index2count:
            vj2index2count[vj] = {i: 1}
        elif i not in vj2index2count[vj]:
            vj2index2count[vj][i] = 1
        else:
            vj2index2count[vj][i] += 1
    
    vj2newclones = {}
    numclone = 0
    for vj, i2count in vj2index2count.iteritems():
        clones = pickle.load(gzip.open(os.path.join(sampledir, vj), 'rb'))
        newclones = []
        for i, count in i2count.iteritems():
            newclone = copy.deepcopy(clones[i])
            newclone.count = count
            newclone.freq = float(count)/size
            newclones.append(newclone)
        vjoutfile = os.path.join(outdir, vj)
        pickle.dump(newclones, gzip.open(vjoutfile, "wb"))
        numclone += len(newclones)
    newsample.size = size
    newsample.numclone = numclone
    outsamfile = os.path.join(outdir, sample.name)
    pickle.dump(newsample, gzip.open(outsamfile, "wb"))

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

#class Sampling(Target):
#    def __init__(self, sample, size, outdir):
#        Target.__init__(self)
#        self.sample = sample
#        self.size = size
#        self.outdir = outdir
#
#    def run(self):
#        subsample = sampling(self.sample, self.size)
#        outfile = os.path.join(self.outdir, "%s.pickle" % sample.name)
#        pickle.dump(subsample, gzip.open(outfile, "wb"))

class SampleAnalysis0(Target):
    '''General child job Obj to do analysis for a specific sample
    '''
    def __init__(self, sample, samdir, outdir, func, *func_args):
        Target.__init__(self)
        self.sample = sample
        self.samdir = samdir
        self.outdir = outdir
        self.func = func
        self.func_args = func_args
    
    def run(self):
        self.func(self.sample, self.samdir, self.outdir, args=self.func_args)

class SampleAnalysis(Target):
    '''General child job Obj to do analysis for a specific sample
    '''
    def __init__(self, sample, samdir, outfile, func, *func_args):
        Target.__init__(self)
        self.sample = sample
        self.samdir = samdir
        self.outfile = outfile
        self.func = func
        self.func_args = func_args
    
    def run(self):
        result_obj = self.func(self.sample, self.samdir, args=self.func_args)
        pickle.dump(result_obj, gzip.open(self.outfile, 'wb')) 

