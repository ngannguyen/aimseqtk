
'''
Make latex table for trackclone results
'''

import sys
import aimseqtk.lib.common as libcommon

def trackclone_latex_tab(f, clone2groups):
    sortedclones = libcommon.sort_by_gene_number(clone2groups.keys())
    #for clone in sorted(clone2groups.keys()):
    for clone in sortedclones:
        cloneitems = clone.split('_')
        v = cloneitems[0]
        seq = cloneitems[1]
        j = cloneitems[2]
        groupfreqs = clone2groups[clone]
        freqstr = ''
        for gfreqs in groupfreqs:
            if not gfreqs:
                freqstr += "& Absent"
            else:
                freqs = [float(freq) for freq in gfreqs.rstrip(',').split(',')]
                sortedfreqs = sorted(freqs, reverse=True)
                sortedfreqs = ["%.4f" % (freq * 100) for freq in sortedfreqs]
                freqstr += " & %s" % ",".join(sortedfreqs)
        f.write("%s & %s & %s %s\\\\\n" % (v, seq, j, freqstr))
        f.write("\\hline\n")

def trackclone_latex(groups, clone2groups, outfile):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["V", "Sequence", "J"] + groups
    colnames = [c.replace('_', "\_") for c in colnames]
    libcommon.tab_header(f, colnames)
    
    trackclone_latex_tab(f, clone2groups)
    caption = ''
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def read_trackclone_txt(file):
    #TRBV6-4_CASSDTLAADSNEQFF_TRBJ2-1       0.076561
    f = open(file, 'r')
    clone2groups = {}
    firstline = f.readline()
    groups = firstline.rstrip('\n').split('\t')[1:]
    for line in f:
        items = line.rstrip('\n').split('\t')
        clone = items[0]
        clone2groups[clone] = items[1:]
    f.close()
    return groups, clone2groups

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    groups, clone2groups = read_trackclone_txt(infile)
    trackclone_latex(groups, clone2groups, outfile)

if __name__ == '__main__':
    main()

