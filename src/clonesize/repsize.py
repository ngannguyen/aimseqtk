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


def sort_samples_by_group(name2sample, group2names, addgroup, group2avr):
    #Sort the samples by group, within each group, sort by size
    sorted_samples = []
    for group in sorted(group2names.keys()):
        gnames = group2names[group]
        gsamples = [name2sample[n] for n in gnames]
        sorted_gsamples = sorted(gsamples, key=lambda gsam: gsam.size, 
                                           reverse=True) 
        sorted_samples.extend(sorted_gsamples)
        if addgroup:
            assert group in group2avr
            groupavr = group2avr[group]
            sorted_samples.append(groupavr)
    return sorted_samples

def get_group_avr(name2sample, group2names):
    group2avr = {}
    for group, names in group2names.iteritems():
        assert len(names) > 0
        avrname = "%s_Avr" % group
        avrsample = libsample.Sample(avrname, group=group)
        for name in names:
            assert name in name2sample
            sample = name2sample[name]
            avrsample.size += sample.size
            avrsample.numclone += len(sample.clones)
        avrsample.size /= len(names)
        avrsample.numclone /= len(names)
        group2avr[group] = avrsample
    return group2avr

def repsize_latex_tab(f, samples, avrnames=[]):
    for s in samples:
        shaded = ""
        if s.name in avrnames:
            shaded = "\\cellcolor[gray]{0.9} "
        name = s.name.replace('_', '-')
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
    for s in samples:
        f.write("%s\t%d\t%d\n" % (s.name, s.numclone, s.size))
    f.close()

def repsize_table(name2sample, outfile, group2avr={}, group2names={}, tex=False):
    # Sort the samples by group and/or size
    if len(group2avr) > 0:
        avrnames = [gavr.name for gavr in group2avr.values()]
        sortedsamples = sort_samples_by_group(name2sample, group2names, 
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

