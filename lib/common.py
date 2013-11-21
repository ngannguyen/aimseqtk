#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions
'''

import os
import sys


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
    f.write("%s \\\\\n" % (" & ".join(colnames)))
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

