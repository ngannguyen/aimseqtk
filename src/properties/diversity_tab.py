
'''
Make latex table for diversity group comparison results
'''

import sys
import aimseqtk.lib.common as libcommon

def ttest_latex_tab(f, lines):
    for line in lines:
        items = line.rstrip('\n').split('\t')
        index = items[0].replace("_", "\_")
        groups = items[1].split('_')
        pval = items[3]
        mean1 = items[4]
        std1 = items[5]
        mean2 = items[6]
        std2 = items[7]
        f.write("\\bf{%s} & %s & %s & %s & %s $\pm$ %s & %s $\pm$ %s\\\\\n" % 
                 (index, groups[0], groups[1], pval, mean1, std1, mean2, std2))
        f.write("\\hline\n")

def ttest_latex(outfile, lines):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Index", "Group 1", "Group 2", "p value", "Mean 1 $\pm$ Std 1", "Mean 2 $\pm$ Std 2"]
    libcommon.tab_header(f, colnames)

    ttest_latex_tab(f, lines)

    caption = ''
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def read_ttest_txt(infile):
    f = open(infile, 'r')
    header = f.readline()
    lines = f.readlines()
    f.close()
    return lines

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    lines = read_ttest_txt(infile)
    ttest_latex(outfile, lines)

if __name__ == '__main__':
    main()
