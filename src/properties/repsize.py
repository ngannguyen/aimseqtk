#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Make summary table of repertoire size:
<sample name>\t<# of reads>\t<# clones>
'''

import os
import sys

import aimseqtk.lib.sample as libsample
import aimseqtk.lib.common as libcommon


def repsize_latex_tab(f, samples, avrnames=[]):
    for (sname, s) in samples:
        shaded = ""
        if s.name in avrnames:
            shaded = "\\cellcolor[gray]{0.9} "
        name = s.name.replace('_', '\_')
        clone = libcommon.pretty_int(s.numclone)
        size = libcommon.pretty_int(s.size)
        f.write("%s%s & %s%s & %s%s \\\\\n" % (shaded, name, shaded, clone,
                                               shaded, size))
        if s.name in avrnames:
            f.write("\\hline\n")
    if not avrnames:
        f.write("\\hline\n")

def repsize_latex(samples, outfile, avrnames=[]):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Sample", "Clones", "Reads"]
    libcommon.tab_header(f, colnames)
    repsize_latex_tab(f, samples, avrnames)
    caption = ("Repertoire size summary. Columns: `Sample': sample name, "
               + "`Clones': number of unique clones, `Reads': number of "
               + "reads. Rows: different samples, where the shaded rows "
               + "show the average of statistics of each group.")
    label = ''
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def repsize_text(samples, outfile):
    f = open(outfile, 'w')
    f.write("Sample\tClones\tReads\n")
    for (sname, s) in samples:
        f.write("%s\t%d\t%d\n" % (s.name, s.numclone, s.size))
    f.close()

def repsize_table(name2sample, outfile, group2avr={}, group2names={}, tex=False):
    # Sort the samples by group and/or size
    if len(group2avr) > 0:
        avrnames = [gavr[0] for gavr in group2avr.values()]
        sortedsamples = libcommon.sort_objs_by_group(name2sample, group2names,
                                                     True, group2avr)
    else:
        avrnames = []
        sortedsamples = sorted(name2sample.values(), key=lambda s: s.size,
                                                     reverse=True)
    # Write to output file
    if tex:
        repsize_latex(sortedsamples, outfile, avrnames)
    else:
        repsize_text(sortedsamples, outfile)

