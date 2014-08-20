#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Get all possible recombination events that result in an amino acid clone
'''

import os
import re
import sys
import gzip
#import marshal as pickle
import cPickle as pickle
from optparse import OptionGroup

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system

import aimseqtk.lib.clone as lclone
import aimseqtk.lib.sample as lsam
import aimseqtk.lib.common as lcommon


class Devent:
    def __init__(self, d5del, d3del, left_nts, cdr3aa_dstart, cdr3aa_dend, right_nts):
        self.d5del = d5del
        self.d3del = d3del
        self.left_nts = left_nts
        self.cdr3aa_dstart = cdr3aa_dstart
        self.cdr3aa_dend = cdr3aa_dend
        self.right_nts = right_nts

def get_all_codons(aaseq):
    codon_lists = []
    for aa in aaseq:
        codon_lists.append(lcommon.aa2codons(aa))
    return codon_lists

def get_all_nts(codon_lists):
    # recursively divide & conquer
    numlist = len(codon_lists)
    if numlist == 1:
        return codon_lists[0]
    half_index = numlist / 2
    left_lists = codon_lists[: half_index]
    right_lists = codon_lists[half_index: ]
    lnts = get_all_nts(left_lists)
    rnts = get_all_nts(right_lists)
    
    combine_list = []
    for l in lnts:
        for r in rnts:
            combine_list.append(l + r)
    return combine_list

def left_max_match(seq1, seq2):
    # find the max matches starting form base0
    matchseq = ''
    l = min(len(seq1), len(seq2))
    for i in xrange(l):
        if seq1[i] == seq2[i]:
            matchseq += seq1[i]
        else:
            break
    return matchseq

def right_max_match(seq1, seq2):
    # find the max matches staring from the end
    matchseq = ''
    l = min(len(seq1), len(seq2))
    endindex = len(seq1) - l - 1
    for i in xrange(1, l + 1):
        if seq1[-1 * i] == seq2[-1 * i]:
            matchseq = seq1[-1 * i] + matchseq
        else:
            break
    return matchseq

def find_min_vdel(v_nt, cdr3_aa):
    # find the minimum number of V nucleotides need to be deleted to
    # result in the CDR3aa
    v_aa = lcommon.nt2aa(v_nt)
    matchseq = left_max_match(v_aa, cdr3_aa)  # max matched aa seq
    len_ntmatch = len(matchseq) * 3
    min_vdel = len(v_nt) - len_ntmatch
    # see if 1 or 2 of the next v nts can result in the next cdr3 aa
    if min_vdel > 0 and len(cdr3_aa) > len(matchseq):
        if min_vdel >= 2:
            v_nts = v_nt[len_ntmatch: len_ntmatch + 2]
        else:
            v_nts = v_nt[len_ntmatch]
        
        got2nts = False
        firstaa = cdr3_aa[len(matchseq)]
        codons = lcommon.aa2codons(firstaa)
        for codon in codons:
            if re.match(v_nts, codon):
                min_vdel -= len(v_nts)
                got2nts = True
                break
        if len(v_nts) == 2 and not got2nts:
            for codon in codons:
                if re.match(v_nts[0], codon):
                    min_vdel -= 1
                    break
    return min_vdel

def find_min_jdel(j_nt, cdr3_aa):
    # find the minimum number of V nucleotides need to be deleted to
    # result in the CDR3aa
    hang_j = len(j_nt) % 3
    hang_j_nts = j_nt[: hang_j]  # the 5' of j nt seq that not part of a codon
    
    j_aa = lcommon.nt2aa(j_nt[hang_j:])
    matchseq = right_max_match(j_aa, cdr3_aa)
    len_ntmatch = len(matchseq) * 3
    min_jdel = len(j_nt) - len_ntmatch
    if min_jdel > 0 and len(cdr3_aa) > len(matchseq):
        if min_jdel >= 2:
            j_nts = j_nt[min_jdel - 2: min_jdel]
        else:
            j_nts = j_nt[min_jdel - 1]
        
        got2nts = False
        firstaa = cdr3_aa[len(cdr3_aa) - len(matchseq) - 1]
        codons = lcommon.aa2codons(firstaa)
        for codon in codons:
            if codon[-1 * len(j_nts): ] == j_nts:
                min_jdel -= len(j_nts)
                got2nts = True
                break
        if len(j_nts) == 2 and not got2nts:
            for codon in codons:
                if codon[2] == j_nts[1]:
                    min_jdel -= 1
                    break
    return min_jdel

def find_dmatches(d_aa, cdr3_aa):
    matches = []
    for mobj in re.finditer(d_aa, cdr3_aa):
        matches.append((mobj.start(), mobj.end()))
    return matches
    
def find_devents(d_nt, cdr3_aa):
    # find all possible (d5del, d3del) such that the resulted nt can be
    # translated to a subsequence of cdr3_aa
    devents = []  
    dlen = len(d_nt)
    for d5del in xrange(dlen + 1):
        for d3del in xrange(dlen - d5del):
            d_leftover_nt = d_nt[d5del: dlen - d3del]
            dlen2 = len(d_leftover_nt)
            if dlen2 == 0:
                event = Devent(d5del, d3del, '', -1, -1, '')
                devents.append(event)
            elif dlen2 < 3:
                event = Devent(d5del, d3del, d_leftover_nt, -1, -1, '')
                devents.append(event)
                #for left in xrange(dlen2 + 1):
                #    left_nts = d_leftover_nt[: left]
                #    right_nts = d_leftover_nt[left: ]
                #    event = Devent(d5del, d3del, left_nts, -1, -1, right_nts)
                #    devents.append(event)
            else:
                for left in [0, 1, 2]:  # each translation frame
                    d_aa = lcommon.nt2aa(d_leftover_nt[left:])
                    # check to see if d_aa match cdr3_aa
                    if d_aa:  # empty sequence
                        matches = find_dmatches(d_aa, cdr3_aa)
                        if len(matches) > 0:
                            right = left + len(d_aa) * 3
                            left_nts = d_leftover_nt[:left]
                            right_nts = d_leftover_nt[right: ]
                            for match in matches:
                                event = Devent(d5del, d3del, left_nts, match[0],
                                               match[1], right_nts)
                                devents.append(event)
                if len(d_leftover_nt) in [3, 4]:
                    event = Devent(d5del, d3del, d_leftover_nt, -1, -1, '')
                    devents.append(event)
    return devents

def get_vdins_events(vdel, v_nt, devent, cdr3_aa, codonlist=False):
    # check if vdel and devent go togther
    # get all possible vdins nt sequences
    v_cdr3_nt = v_nt if vdel == 0 else v_nt[: -1 * vdel]
    lastvpos = len(v_cdr3_nt) - 1

    if devent.cdr3aa_dstart > 0:  # there is still part of D
        d5pos = devent.cdr3aa_dstart * 3 - len(devent.left_nts)
        if d5pos <= lastvpos:
            return None
        else:
            v_hang = len(v_cdr3_nt) % 3
            vd_ins_len = d5pos - lastvpos - 1  # number of nts between lastvpos and firstd5pos
            vd_inframe = v_hang + vd_ins_len + len(devent.left_nts)
            if vd_inframe % 3 > 0:
                return None
            else:
                vd_inframe_aa = cdr3_aa[len(v_cdr3_nt) / 3: devent.cdr3aa_dstart]
                if not vd_inframe_aa:  # no inserted nts
                    if codonlist:
                        return [[v_cdr3_nt[-1]]]
                    else:
                        return [v_cdr3_nt[-1]]
                else:
                    codon_lists = get_all_codons(vd_inframe_aa)
                    if devent.left_nts:
                        len_left_nts = len(devent.left_nts)
                        d_codons = []
                        for codon in codon_lists[-1]:
                            if codon[-1 * len_left_nts: ] == devent.left_nts:
                                d_codons.append(codon[: -1 * len_left_nts])  # remove d_hang
                        if not d_codons:
                            return None
                        codon_lists[-1] = d_codons
                    if v_hang > 0:
                        v_hang_nts = v_cdr3_nt[-1 * v_hang: ]
                        v_codons = []
                        for codon in codon_lists[0]:
                            if re.match(v_hang_nts, codon):
                                if v_hang == 1:
                                    v_codons.append(codon)  # leave the last v nt
                                else:
                                    v_codons.append(codon[1:])
                        if not v_codons:
                            return None
                        codon_lists[0] = v_codons
                    else:  # add the last v nt
                        last_v_nt = v_cdr3_nt[-1]
                        codon_lists = [[last_v_nt]] + codon_lists

                    if codonlist:
                        return codon_lists
                    else:
                        vdins_nts = get_all_nts(codon_lists)
                        return vdins_nts
    else:
        assert devent.cdr3aa_dend == -1
        return None

def get_djins_events(jdel, j_nt, devent, cdr3_aa, codonlist=False):
    # check if jdel and devent go togther
    # get all possible djins nt sequences
    #print "jdel: %d" % jdel
    j_cdr3_nt = j_nt if jdel == 0 else j_nt[jdel: ]
    #print "j_cdr3_nt: %s" % j_cdr3_nt
    firstjpos = len(cdr3_aa) * 3 - len(j_cdr3_nt)

    if devent.cdr3aa_dend > 0:  # there is still part of D
        d3pos = devent.cdr3aa_dend * 3 + len(devent.right_nts) - 1  # inclusive
        if d3pos >= firstjpos:  # overlap d & j
            return None
        else:
            j_hang = len(j_cdr3_nt) % 3
            #print j_hang  # 2

            dj_ins_len = firstjpos - d3pos - 1  # number of nts between lastvpos and firstd5pos
            dj_inframe = len(devent.right_nts) + dj_ins_len + j_hang
            if dj_inframe % 3 > 0:
                return None
            else:
                dj_inframe_aa = cdr3_aa[devent.cdr3aa_dend: len(cdr3_aa) - len(j_cdr3_nt) / 3]
                #print dj_inframe_aa  # VN
                if not dj_inframe_aa:  # no inserted nts
                    if codonlist:
                        return [[j_cdr3_nt[0]]]
                    else:
                        return [j_cdr3_nt[0]]
                else:
                    codon_lists = get_all_codons(dj_inframe_aa)
                    print codon_lists
                    if devent.right_nts:
                        len_right_nts = len(devent.right_nts)
                        d_codons = []
                        for codon in codon_lists[0]:
                            if codon[: len_right_nts] == devent.right_nts:
                                d_codons.append(codon[len_right_nts: ])  # remove d_hang
                        if not d_codons:
                            return None
                        codon_lists[0] = d_codons
                    if j_hang > 0:
                        j_hang_nts = j_cdr3_nt[: j_hang]
                        j_codons = []
                        for codon in codon_lists[-1]:
                            if codon[-1 * j_hang: ] == j_hang_nts:
                                if j_hang == 1:
                                    j_codons.append(codon)  # leave the last v nt
                                else:
                                    j_codons.append(codon[: -1])
                        if not j_codons:
                            return None
                        codon_lists[-1] = j_codons
                    else:  # add the last v nt
                        first_j_nt = j_cdr3_nt[0]
                        codon_lists = codon_lists + [[first_j_nt]]

                    if codonlist:
                        return codon_lists
                    else:
                        djins_nts = get_all_nts(codon_lists)
                        return djins_nts
    else:
        assert devent.cdr3aa_dstart == -1
        return None

def get_vjins_emptyd(v_nt, vdel, j_nt, jdel, d_nts, cdr3_aa):
    v_cdr3_nt = v_nt if vdel == 0 else v_nt[: -1 * vdel]
    lastvpos = len(v_cdr3_nt) - 1
    j_cdr3_nt = j_nt if jdel == 0 else j_nt[jdel: ]
    firstjpos = len(cdr3_aa) * 3 - len(j_cdr3_nt)
    
    if lastvpos >= firstjpos:  # overlap
        return None

    v_hang = len(v_cdr3_nt) % 3
    j_hang = len(j_cdr3_nt) % 3
    vj_inframe_aa = cdr3_aa[len(v_cdr3_nt) / 3: len(cdr3_aa) - len(j_cdr3_nt) / 3] 
    if not vj_inframe_aa:  # no insertion base
        if not d_nts:
            return [v_cdr3_nt[-1] + j_cdr3_nt[0]]  #last_v and first_j
        else:
            return None

    codon_lists = get_all_codons(vj_inframe_aa)
    if v_hang > 0:
        v_hang_nts = v_cdr3_nt[-1 * v_hang: ]
        v_codons = []
        for codon in codon_lists[0]:
            if re.match(v_hang_nts, codon):
                if v_hang == 1:
                    v_codons.append(codon)  # leave the last v nt
                else:
                    v_codons.append(codon[1:])
        assert v_codons
        codon_lists[0] = v_codons
    else:
        last_v_nt = v_cdr3_nt[-1]
        codon_lists = [[last_v_nt]] + codon_lists

    if j_hang > 0:
        j_hang_nts = j_cdr3_nt[: j_hang]
        j_codons = []
        for codon in codon_lists[-1]:
            if codon[-1 * j_hang: ] == j_hang_nts:
                if j_hang == 1:
                    j_codons.append(codon)  # leave the last v nt
                else:
                    j_codons.append(codon[: -1])
        if not j_codons:
            return None
        codon_lists[-1] = j_codons
    else:  # add first j nt
        first_j_nt = j_cdr3_nt[0]
        codon_lists = codon_lists + [[first_j_nt]]

    vjins_nts = get_all_nts(codon_lists)
    vjins_has_d = []
    for vjins in vjins_nts:
        if re.search(d_nts, vjins):
            vjins_has_d.append(vjins)
    return vjins_has_d

class Get_Vjins(Target):
    def __init__(self, clone, v_nt, min_vdel, max_vdel, j_nt, min_jdel,
                       max_jdel, d, devent, cdr3_aa, outfile):
        Target.__init__(self)
        self.clone = clone
        self.v_nt = v_nt
        self.min_vdel = min_vdel
        self.max_vdel = max_vdel
        self.j_nt = j_nt
        self.min_jdel = min_jdel
        self.max_jdel = max_jdel
        self.d = d
        self.devent = devent
        self.cdr3_aa = cdr3_aa
        self.outfile = outfile

    def run(self):
        items = self.clone.split('_')
        v = items[0]
        j = items[2]
        events = []
        batchsize = 100000
        currbatch = 0

        for vdel in xrange(self.min_vdel, self.max_vdel + 1):
            v_cdr3_nt = self.v_nt if vdel == 0 else self.v_nt[: -1 * vdel]
            v_hang = len(v_cdr3_nt) % 3
            for jdel in xrange(self.min_jdel, self.max_jdel + 1):
                j_cdr3_nt = self.j_nt if jdel == 0 else self.j_nt[jdel: ]
                d_nts = self.devent.left_nts + self.devent.right_nts
                vjins_nts = get_vjins_emptyd(self.v_nt, vdel, self.j_nt, jdel,
                                             d_nts, self.cdr3_aa)
                if vjins_nts is None:
                    continue
                
                #self.logToMaster("Empty D: vdel: %d, jdel: %d, vjins: %d\n" % (vdel, jdel, len(vjins_nts)))

                for vjins in vjins_nts:
                    assert len(vjins) >= 2  # because it has the last v and the first j
                    for mobj in re.finditer(d_nts, vjins):
                        start = mobj.start()
                        end = mobj.end()
                        # check special case
                        v_hang_right = 0 if v_hang == 0 else 3 - v_hang
                        if len(d_nts) == 3:
                            if (start - v_hang_right - 1) % 3 == 0:
                                continue
                        elif len(d_nts) == 4:
                            if (start - v_hang_right - 1) != 1:
                                continue
                        
                        vdins = vjins[: start]
                        djins = vjins[end: ]
                        #assert len(vdins) + len(djins) + len(d_nts) == len(vjins)
                        cdr3_nt = v_cdr3_nt + vjins[1:-1] + j_cdr3_nt
                        assert lcommon.nt2aa(cdr3_nt) == self.cdr3_aa 
                        event = lclone.Cdr3Clone(1, cdr3_nt, v, j, d=self.d,
                              aa=self.cdr3_aa, vdel=vdel, jdel=jdel,
                              d5del=self.devent.d5del, d3del=self.devent.d3del,
                              vdins=vdins, djins=djins)
                        events.append(event)

                        if len(events) >= batchsize:
                            outfile = "%s_%d" % (self.outfile, currbatch)
                            pickle.dump(events, gzip.open(outfile, 'wb'))
                            currbatch += 1
                            events = []
        if len(events) > 0:
            outfile = "%s_%d" % (self.outfile, currbatch)
            pickle.dump(events, gzip.open(outfile, 'wb'))

class Get_Ins(Target):
    def __init__(self, func, ntdel, nt, devent, cdr3_aa, outfile):
        Target.__init__(self)
        self.func = func
        self.ntdel = ntdel
        self.nt = nt
        self.devent = devent
        self.cdr3_aa = cdr3_aa
        self.outfile = outfile

    def run(self):
        ins_nts = self.func(self.ntdel, self.nt, self.devent, self.cdr3_aa)
        if ins_nts is not None:
            self.logToMaster("get_ins: %d" % len(ins_nts))
            pickle.dump(ins_nts, gzip.open(self.outfile, 'wb'))

class Get_Vd_Dj_Ins_Agg(Target):
    def __init__(self, clone, vdir, jdir, v_nt, j_nt, d, d_nt, devent, outfile):
        Target.__init__(self)
        self.clone = clone
        self.vdir = vdir
        self.jdir = jdir
        self.v_nt = v_nt
        self.j_nt = j_nt
        self.d = d
        self.d_nt = d_nt
        self.devent = devent
        self.outfile = outfile

    def run(self):
        events = []
        batchsize = 100000
        currbatch = 0
        items = self.clone.split('_')
        v = items[0]
        cdr3_aa = items[1]
        j = items[2]
        if self.devent.d3del == 0:
            d_cdr3_nt = self.d_nt[self.devent.d5del: ]
        else:
            d_cdr3_nt = self.d_nt[self.devent.d5del: -1 * self.devent.d3del]

        for vdel in os.listdir(self.vdir):
            vfile = os.path.join(self.vdir, vdel)
            vd_ins_nts = pickle.load(gzip.open(vfile, 'rb'))
            vdel = int(vdel)
            self.logToMaster("vdel: %d; vd_ins_nts: %d" % (vdel, len(vd_ins_nts)))
            v_cdr3_nt = self.v_nt if vdel == 0 else self.v_nt[: -1 * vdel]
            for jdel in os.listdir(self.jdir):
                jfile = os.path.join(self.jdir, jdel)
                dj_ins_nts = pickle.load(gzip.open(jfile, 'rb'))
                jdel = int(jdel)
                self.logToMaster("jdel: %d; dj_ins_nts: %d" % (jdel, len(dj_ins_nts)))
                j_cdr3_nt = self.j_nt if jdel == 0 else self.j_nt[jdel: ]
                
                for vd_ins in vd_ins_nts:
                    for dj_ins in dj_ins_nts:
                        cdr3_nt = (v_cdr3_nt + vd_ins[1: ] + d_cdr3_nt +
                                   dj_ins[: -1] + j_cdr3_nt)
                        if lcommon.nt2aa(cdr3_nt) != cdr3_aa:
                            print cdr3_nt
                            print lcommon.nt2aa(cdr3_nt)
                            print cdr3_aa
                        assert lcommon.nt2aa(cdr3_nt) == cdr3_aa
                        event = lclone.Cdr3Clone(1, cdr3_nt, v, j, d=self.d,
                              aa=cdr3_aa, vdel=vdel, jdel=jdel,
                              d5del=self.devent.d5del, d3del=self.devent.d3del,
                              vdins=vd_ins, djins = dj_ins)
                        events.append(event)
                        if len(events) >= batchsize:
                            outfile = "%s_%d" % (self.outfile, currbatch)
                            pickle.dump(events, gzip.open(outfile, 'wb'))
                            currbatch += 1
                            events = []
        if len(events) > 0:
            outfile = "%s_%d" % (self.outfile, currbatch)
            pickle.dump(events, gzip.open(outfile, 'wb'))

class Get_Vd_Dj_Ins(Target):
    def __init__(self, clone, v_nt, min_vdel, max_vdel, j_nt, min_jdel,
                       max_jdel, d, d_nt, devent, cdr3_aa, outdir):
        Target.__init__(self)
        self.clone = clone
        self.v_nt = v_nt
        self.min_vdel = min_vdel
        self.max_vdel = max_vdel
        self.j_nt = j_nt
        self.min_jdel = min_jdel
        self.max_jdel = max_jdel
        self.d = d
        self.d_nt = d_nt
        self.devent = devent
        self.cdr3_aa = cdr3_aa
        self.outdir = outdir

    def run(self):
        vdir = os.path.join(self.outdir, "vdels")
        system("mkdir -p %s" % vdir)
        for vdel in xrange(self.min_vdel, self.max_vdel + 1):
            voutfile = os.path.join(vdir, str(vdel))
            self.addChildTarget(Get_Ins(get_vdins_events, vdel,
                        self.v_nt, self.devent, self.cdr3_aa, voutfile))
        jdir = os.path.join(self.outdir, 'jdels')
        system("mkdir -p %s" % jdir)
        for jdel in xrange(self.min_jdel, self.max_jdel + 1):
            joutfile = os.path.join(jdir, str(jdel))
            self.addChildTarget(Get_Ins(get_djins_events, jdel,
                        self.j_nt, self.devent, self.cdr3_aa, joutfile))
        
        outfile = os.path.join(self.outdir, "events")
        self.setFollowOnTarget(Get_Vd_Dj_Ins_Agg(self.clone, vdir, jdir,
                self.v_nt, self.j_nt, self.d, self.d_nt, self.devent, outfile))

class GetCloneEvents(Target):
    def __init__(self, clone, aaseq, vseq, jseq, d2seq, outdir):
        Target.__init__(self)
        self.clone = clone
        self.aaseq = aaseq
        self.vseq = vseq
        self.jseq = jseq
        self.d2seq = d2seq
        self.outdir = outdir
    
    def run(self):
        self.logToMaster("Getting recomb. events for clone %s ..." % self.clone)
        max_vdel = len(self.vseq) - 3
        min_vdel = find_min_vdel(self.vseq, self.aaseq)
        max_jdel = len(self.jseq) - 3
        min_jdel = find_min_jdel(self.jseq, self.aaseq)
        self.logToMaster("Vdel: <%d-%d>" % (min_vdel, max_vdel))
        self.logToMaster("Jdel: <%d-%d>" % (min_jdel, max_jdel))

        for d, dseq in self.d2seq.iteritems():
            devents = find_devents(dseq, self.aaseq)
            self.logToMaster("%d number of devents" % (len(devents)))
            # DEBUG
            #numempty = 0
            #for devent in devents:
            #    if devent.cdr3aa_dstart == -1:
            #        numempty += 1
            #self.logToMaster("\t%d empty D, %d non_empty_D\n" % (numempty, len(devents) - numempty))
            # END DEBUG
            for i, devent in enumerate(devents):
                outdir = os.path.join(self.outdir, d, str(i))  #outdir/clone/d/i
                system("mkdir -p %s" % outdir)
                
                if devent.cdr3aa_dstart == -1:
                    dempty_file = os.path.join(outdir, "d_empty")
                    self.addChildTarget(Get_Vjins(self.clone, self.vseq,
                                min_vdel, max_vdel, self.jseq, min_jdel,
                                max_jdel, d, devent, self.aaseq, dempty_file))
                else:
                    self.addChildTarget(Get_Vd_Dj_Ins(self.clone, self.vseq,
                                min_vdel, max_vdel, self.jseq, min_jdel,
                                max_jdel, d, dseq, devent, self.aaseq, outdir))
        self.setFollowOnTarget(CloneEventsAgg(self.outdir))

class CloneEventsAgg(Target):
    def __init__(self, outdir):
        Target.__init__(self)
        self.outdir = outdir

    def run(self):
        clone2numevents = {}
        for c in os.listdir(self.outdir):
            clonedir = os.path.join(self.outdir, c)
            count = 0
            for d in os.listdir(clonedir):
                ddir = os.path.join(clonedir, d)
                for devent in os.listdir(ddir):
                    devent_dir = os.path.join(ddir, devent)
                    for file in os.listdir(devent_dir):
                        filepath = os.path.join(devent_dir, file)
                        if not os.path.isdir(filepath):
                            events = pickle.load(gzip.open(filepath, 'rb'))
                            count += len(events)
            clone2numevents[c] = count

        outfile = os.path.join(self.outdir, "summary.txt")
        f = open(outfile, 'w')
        f.write("#Clone\tNumber of recombination events\n")
        for c, count in clone2numevents.iteritems():
            f.write("%s\t%d\n" % (c, count))

class Setup(Target):
    def __init__(self, args):
        Target.__init__(self)
        self.clonefile = args[0]
        self.vfile = args[1]
        self.jfile = args[2]
        self.dfile = args[3]
        self.outdir = args[4]

    def run(self):
        self.logToMaster("Setting up...\n")
        if not os.path.exists(self.outdir):
            system("mkdir -p %s" % self.outdir)

        #global_dir = self.getGlobalTempDir()
        clones = lcommon.read_list(self.clonefile)
        v2seq = lcommon.read_dict(self.vfile, cap=True)
        j2seq = lcommon.read_dict(self.jfile, cap=True)
        d2seq = lcommon.read_dict(self.dfile, cap=True)
        self.logToMaster("Done processing the input files.\n")

        for clone in clones:
            items = clone.split('_')
            v = items[0]
            seq = items[1]
            j = items[2]

            vseq = v2seq[v]
            jseq = j2seq[j]
            outdir = os.path.join(self.outdir, clone)
            system("mkdir -p %s" % outdir)
            self.addChildTarget(GetCloneEvents(clone, seq, vseq, jseq, d2seq,
                                               outdir))

def main():
    usage = ("%prog <clonefile> <vfile> <jfile> <dfile> <outdir>")
    parser = lcommon.init_options(usage)
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    i = Stack(Setup(args)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from aimseqtk.src.recomb.aa_events import *
    main()
