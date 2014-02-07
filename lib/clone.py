#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Object represents a unique clonotype and related functions
'''


class Cdr3Clone():
    def __init__(self, count, nuc, v, j, d='', aa=None, sample=None,
                 normcount=None, freq=None, vdel=None, jdel=None, d5del=None,
                 d3del=None, vdins='', djins='', id=None):
        self.count = count
        self.nuc = nuc
        self.v = v
        self.j = j
        self.d = d
        self.vdel = vdel
        self.jdel = jdel
        self.d5del = d5del
        self.d3del = d3del
        self.vdins = ''
        self.djins = ''
        self.id = id
        self.sample = sample
        self.freq = freq
        self.aa = aa
        if normcount:
            self.normcount = normcount
        else:
            self.normcount = count
         
    def __getitem__(self, attr):
        if attr not in self.__dict__:
            raise KeyError("Cdr3Clone does not have attribute %s" % attr)
        return self.__dict__[attr]

    def __setitem__(self, attr, val):
        self.__dict__[attr] = val

    def getitems(self):
        return self.__dict__.keys()

    def get_vseqj(self, nuc=False):
        if nuc:
            seq = self.nuc
        else:
            seq = self.aa
        return "%s_%s_%s" % (self.v, seq, self.j)

class Clone():
    def __init__(self, count, freq, nuc, vgenes, jgenes, dgenes=[],
                 valleles=[], dalleles=[], jalleles=[], 
                 cdr3nuc=None, cdr3aa=None, aa=None, productive=None, 
                 lastvpos=None, firstdpos=None, lastdpos=None, firstjpos=None, 
                 vdel=None, d5del=None, d3del=None, jdel=None, id=None):
        
        # Basic information
        self.count = count  # Read count
        self.normcount = count
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
        self.firstdpos = firstdpos  # first position of D segment
        self.lastdpos = lastdpos  # last position of D segment
        self.firstjpos = firstjpos  # first position of J segment
        self.vdel = vdel  # number of bases of the vsegment got deleted
        self.d5del = d5del
        self.d3del = d3del 
        self.jdel = jdel
        self.id = id  # clone id
        self.vjseq_ids = []
        
    def __getitem__(self, name):
        if name not in self.__dict__:
            raise KeyError("Clone does not have attribute %s" % name)
        return self.__dict__[name]

    def __setitem__(self, name, val):
        self.__dict__[name] = val

    def getitems(self):
        return self.__dict__.keys()

    def get_sorted_items(self):
        items = ['count', 'freq', 'nuc', 'vgenes', 'jgenes', 'dgenes',
                 'cdr3nuc', 'valleles', 'dalleles', 'jalleles', 'cdr3aa', 'aa',
                 'productive', 'lastvpos', 'firstdpos', 'lastdpos',
                 'firstjpos', 'vdel', 'd5del', 'd3del', 'jdel', 'id']
        return items
    
    def getstr(self):
        fields = self.get_sorted_items()
        vals = []
        present_fields = self.getitems()
        for field in fields:
            if field not in present_fields or self[field] is None:
                vals.append('')
            else:
                val = self[field]
                if isinstance(val, list):
                    vals.append(", ".join(val))
                else:
                    vals.append(str(val))
        return "\t".join(vals)

    def set_normcount(self, normcount):
        #self.normcount = normcount
        self.__dict__['normcount'] = normcount

    def get_vjseq_ids(self):
        ids = []
        seq = self.cdr3nuc
        if self.cdr3nuc is None:
            seq = self.nuc
        for v in self.vgenes:
            for j in self.jgenes:
                id = "%s_%s_%s" % (v, seq, j)
                ids.append(id)
        self.vjseq_ids = ids
        return ids

#======= Read and write functions ==========
def clone_columns():
    cols = ['count', 'freq', 'nuc', 'vgenes', 'jgenes', 'dgenes',
             'cdr3nuc', 'valleles', 'dalleles', 'jalleles', 'cdr3aa', 'aa',
             'productive', 'lastvpos', 'firstdpos', 'lastdpos',
             'firstjpos', 'vdel', 'd5del', 'd3del', 'jdel', 'id']
    return cols

def clone_parseline(line, index2col):
    items = line.strip().split('\t')
    if len(items) != len(index2col):
        sys.stderr.write(("Inconsistent number of columns between the " +
                          "following line and the header line, skipped it:" +
                          "\nLine:\n%s\n" % line))
        return None
    col2val = {}
    valid_cols = clone_columns()
    for i, col in index2col.iteritems():
        if col in valid_cols:
            col2val[col] = items[i]
    
    # Return None if line does not have minimum required fields.
    required_cols = ['count', 'freq', 'nuc', 'vgenes', 'jgenes']
    for c in required_cols:
        if c not in col2val or not col2val[c]:
            return None
    count = int(col2val['count'])
    freq = float(col2val['freq'])
    nuc = col2val['nuc']
    vgenes = col2val['vgenes'].split(', ')
    jgenes = col2val['jgenes'].split(', ')
    
    clone = Clone(count, freq, nuc, vgenes, jgenes)

    list_cols = ['dgenes', 'valleles', 'dalleles', 'jalleles']
    for c in list_cols:
        if c in col2val and col2val[c]:
            clone[c] = col2val[c].split(', ')
    str_cols = ['cdr3nuc', 'cdr3aa', 'aa', 'id']
    for c in str_cols:
        if c in col2val and col2val[c]:
            clone[c] = col2val[c]
    int_cols = ['lastvpos', 'firstdpos', 'lastdpos', 'firstjpos',
                'vdel', 'd5del', 'd3del', 'jdel']
    for c in int_cols:
        if c in col2val and col2val[c]:
            clone[c] = int(col2val[c])
    if 'productive' in col2val and col2val['productive']:
        if col2val['productive'].lower() == 'true':
            clone.productive = True
        else:
            clone.productive = False
    return clone
