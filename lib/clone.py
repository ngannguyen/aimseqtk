#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a unique clonotype and related functions
'''

class Clone():
    def __init__(self, count, freq, nuc, vgenes, jgenes, dgenes=[],
                 valleles=[], dalleles=[], jalleles=[], 
                 cdr3nuc=None, cdr3aa=None, aa=None, productive=None, 
                 lastvpos=None, firstdpos=None, lastdpos=None, firstjpos=None, 
                 vdel=None, d5del=None, d3del=None, jdel=None, id=None):
        
        # Basic information
        self.count = count  # Read count
        self.freq = freq  # Read count/Total read count
        self.nuc = nuc  # full nt seq (cdr3nuc is a subseq of nuc)
        self.vgenes = vgenes  # list of v genes, e.g [TRBV2.7]
        self.jgenes = jgenes  # list of j genes
        self.dgenes = dgenes  # list of d genes
        
        # Optional information
        self.cdr3nuc = cdr3nuc  # cdr3 nucleotide sequence
        self.valleles = valleles  # e.g: [TRBV2.7*01]
        self.dalleles = dalleles
        self.jalleles = jalleles
        self.cdr3aa = cdr3aa  # cdr3 aa seq, commonly starts wt C and ends wt F
        self.aa = aa  # full aa seq
        self.productive = productive  # True if the clone is productive
        
        # lastvpos, firstdpos, lastdpos and firstjpos are relative to "nuc"
        # all are inclusive, base 0
        self.lastvpos = lastvpos  # last position of V segment
        self.fistdpos = firstdpos  # first position of D segment
        self.lastdpos = lastdpos  # last position of D segment
        self.firstjpos = firstjpos  # first position of J segment
        self.vdel = vdel  # number of bases of the vsegment got deleted
        self.d5del = d5del
        self.d3del = d3del 
        self.jdel = jdel
        self.id = id  # clone id
        
    

