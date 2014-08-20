
'''
Make latex table for similarity group comparison results
'''

import sys
import aimseqtk.lib.common as libcommon

def group_cmp_latex_tab(f, lines):
    for line in lines:
        items = line.rstrip('\n').split('\t')
        groups = items[0].replace("_", "\_").split(', ')
        pval = items[1]
        mean1 = items[3]
        mean2 = items[4]
        f.write("%s & %s & %s & %s & %s\\\\\n" % 
                 (groups[0], groups[1], pval, mean1, mean2))
        f.write("\\hline\n")

def group_cmp_latex(outfile, lines):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Category 1", "Category 2", "p value", "Mean 1 $\pm$ Std 1", "Mean 2 $\pm$ Std 2"]
    libcommon.tab_header(f, colnames)

    group_cmp_latex_tab(f, lines)

    caption = ''
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def read_group_cmp_txt(infile):
    f = open(infile, 'r')
    header = f.readline()
    lines = f.readlines()
    f.close()
    return lines

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    lines = read_group_cmp_txt(infile)
    group_cmp_latex(outfile, lines)

if __name__ == '__main__':
    main()
