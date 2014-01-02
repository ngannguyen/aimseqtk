#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys
import copy
import numbers
import cPickle as pickle
from optparse import OptionParser

from jobTree.scriptTree.target import Target


class InputError(Exception):
    pass

def get_cumulative(vals, forward=False):
    if forward:
        return [sum(vals[: i + 1]) for i in xrange(len(vals))]
    else:
        return [sum(vals[i: ]) for i in xrange(len(vals))]

def get_pc(myval, total):
    if total == 0:
        return 0.0
    else:
        return 100.0 * myval / total

def sort_dict_by_value(mydict, keyfunc=None):
    # Return a list of tuples
    items = [(v, k) for k, v in mydict.items()]
    if keyfunc is None:
        items = sorted(items, reverse=True)
    else:
        items = sorted(items, reverse=True, key=keyfunc)
    return [(k, v) for v, k in items] 

def sort_objs_by_group(name2obj, group2names, addgroup, group2avr,
                                                        keyfunc=None):
    # Sort the objs by group, within each group, sort by keyfunc if
    # specified, else by default
    sorted_items = []
    for group in sorted(group2names.keys()):
        gnames = group2names[group]
        gname2obj = dict((n, name2obj[n]) for n in gnames)
        gsorted_items = sort_dict_by_value(gname2obj, keyfunc=keyfunc)
        sorted_items.extend(gsorted_items)
        if addgroup:
            assert group in group2avr
            groupavr = group2avr[group]
            sorted_items.append(groupavr)
    return sorted_items

def get_group_avr(name2obj, group2names):
    # for each group, return an ojb whose attributes are average of the
    # group member objs attributes. Note: obj must have func "getitems",
    # which return the obj attributes, defined
    group2avr = {}
    if not group2names:
        return group2avr

    for group, names in group2names.iteritems():
        assert len(names) > 0
        avrname = "%s_Avr" % group
        avrobj = copy.copy(name2obj[names[0]]) 
        # Reset the avrobj:
        for attr in avrobj.getitems():
            val = avrobj[attr]
            if isinstance(val, numbers.Number):
                avrobj.__setitem__(attr, 0)
            elif isinstance(val, list):
                if len(val) > 0:
                    numericlist = True
                    for vitem in val:
                        if not isinstance(vitem, numbers.Number):
                            numericlist = False
                            break
                    if numericlist:
                        init_val = [0 for vitem in val]
                        avrobj.__setitem__(attr, init_val)
                    else:
                        avrobj.__setitem__(attr, None)
            elif attr == 'name':
                avrobj.__setitem__(attr, avrname)
            elif attr == 'group':
                avrobj.__setitem__(attr, group)
            elif not isinstance(val, str):
                avrobj.__setitem__(attr, None)

        # Update all the numeric attributes
        for name in names:  # each member
            assert name in name2obj
            obj = name2obj[name]
            for attr in obj.getitems():
                val = obj[attr]
                if isinstance(val, numbers.Number):
                    newval = avrobj[attr] + val
                    avrobj.__setitem__(attr, newval)
                elif isinstance(val, list):  # attr is a list of numeric values
                    if len(val) > 0:
                        for i, vitem in enumerate(val):
                            if isinstance(vitem, numbers.Number):
                                newvitem = avrobj[attr][i] + vitem
                                avrobj[attr][i] = newvitem

        for attr in avrobj.getitems():
            val = avrobj[attr]
            if isinstance(val, numbers.Number):
                avrval = val/len(names)
                avrobj.__setitem__(attr, avrval)
            elif isinstance(val, list):
                if len(val) > 0:
                    for i, vitem in enumerate(val):
                        if isinstance(vitem, numbers.Number):
                            avrvitem = vitem/len(names)
                            avrobj[attr][i] = avrvitem

        group2avr[group] = (avrname, avrobj)
    return group2avr

def load_pickledir(pickledir):
    objs = []
    for file in os.listdir(pickledir):
        filepath = os.path.join(pickledir, file)
        obj = pickle.load(gzip.open(filepath, "rb"))
        objs.append(obj)
    return objs

def load_pickledir_to_dict(pickledir):
    name2obj = {}
    for file in os.listdir(pickledir):
        filepath = os.path.join(pickledir, file)
        name = file.splitext()[0]
        obj = pickle.load(gzip.open(filepath, "rb"))
        name2obj[name] = obj
    return name2obj

def get_val2keys(key2vals):
    val2keys = {}
    for key, vals in key2vals.iteritems():
        for val in vals:
            if val not in val2keys:
                val2keys[val] = [key]
            else:
                val2keys[val].append(key)
    return val2keys

def get_val2key_1to1(key2vals):
    # revert a dictionary, requires that 1 value corresponds to only 1 key
    val2key = {}
    for k, vals in key2vals.iteritems():
        for v in vals:
            if v in val2key:
                raise InputError(("%s belongs to multiple groups: " % v
                                  + "%s, %s" % (k, val2key[v]))) 
            val2key[v] = k
    return val2key

def read_list(file):
    items = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) > 0 and line[0] != '#':
            items.append(line)
    f.close()
    return items

def check_options_dir(dir):
    if not os.path.exists(dir):
        raise InputError("Directory %s does not exist." % dir)
    if not os.path.isdir(dir):
        raise InputError("%s is not a directory." % dir)

def check_options_file(file):
    if not os.path.exists(file):
        raise InputError("File %s does not exists." % file)

def init_options(usage=None):
    if usage is None:
        usage = "%prog [options]"
    parser = OptionParser(usage = usage)
    return parser

def get_index2item(line):
    items = line.rstrip("\n").split("\t")
    index2item = {}
    for i, item in enumerate(items):
        index2item[i] = item
    return index2item

def nt2aa(nt):
    # Translate nucleotide sequence to amino acid sequence
    codon2aa = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': 'Z', 'TGA': 'Z',
    'TTG': 'L', 'TCG': 'S', 'TAG': 'Z', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' 
    }
    aaseq = ''
    for i in xrange(0, len(nt)/3):
        codon = nt[i*3: i*3+3].upper()
        aaseq = aaseq + codon2aa[codon]
    return aaseq

def rc(nt):
    #Return a reverse complement of the input sequence.
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    rcnt = ''
    for i in xrange( len(nt) -1, -1, -1):
        if nt[i] in complement:
            rcnt += complement[nt[i]] 
        else:
            rcnt += nt[i]
    return rcnt

def soft_float(mystr):
    # If possible, convert mystr into a float, otherwise return mystr
    try:
        myfloat = float(mystr)
    except:
        myfloat = mystr
    return myfloat

def soft_int(mystr):
    # If possible, convert mystr into an int, otherwise return mystr
    try:
        myfloat = float(mystr)
        myint = int(myfloat)
    except:
        myint = mystr
    return myint

#================ COMMON STRUCTUREs ==========
class Analysis(Target):
    '''Parent object of each type of analyses, e.g Geneusage, Lendist,
    Similarity, Diversity etc
    '''
    def __init__(self, samples, outdir, opts=None, *args):
        Target.__init__(self)
        self.samples = samples
        self.outdir = outdir
        self.opts = opts
        self.args = args

class StatAnalyses(Target):
    '''Parent object which reads in stat objs (one/sample) from indir
    '''
    def __init__(self, indir, outdir, opts=None, *args):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.opts = opts
        self.args = args
        self.name2obj = load_pickledir_to_dict(indir)

#============ LATEX related function ============
def pretty_int(number):
    numstr = str(number)
    prettyStr = ''
    for i in xrange(len(numstr)):
        if i > 0 and (len(numstr) - i) % 3 == 0:
            prettyStr = prettyStr + ','
        prettyStr = prettyStr + numstr[i]
    return prettyStr

def pretty_float(num):
    if num == 0:
        return "0"
    else:
        return "%.2e" % num

def tab_header(f, colnames):
    f.write("\\begin{table}\n")
    f.write("\\centering\n")
    f.write("\\scalebox{1.0}{%\n")
    f.write("\\begin{tabular}{c%s}\n" % ("|c" * (len(colnames) - 1)))
    f.write("\\hline\n")
    boldedcols = ["\\bf{%s}" % c for c in colnames]
    f.write("%s \\\\\n" % (" & ".join(boldedcols)))
    f.write("\\hline\n")

def write_doc_start(f):
    f.write("\\documentclass[11pt]{article}\n")
    f.write("\\usepackage{epsfig}\n")
    f.write("\\usepackage{multirow}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{array}\n")
    f.write("\\usepackage{color, colortbl}\n")
    f.write("\\usepackage{rotating}\n")
    f.write("\\usepackage[table]{xcolor}\n")
    f.write("\\definecolor{LightCyan}{rgb}{0.88,1,1}")
    f.write("\n" )

    f.write("\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n")
    f.write("\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n")
    f.write("\n")

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")

    f.write("\\begin{document}\n")
    return

def write_doc_end(f):
    f.write( "\\end{document}\n" )

def tab_closer(f, caption, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %caption)
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def sideways_tab_closer(f, caption, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %caption)
    f.write("\\label{%s}" %label)
    f.write("\\end{sidewaystable}\n\n")

