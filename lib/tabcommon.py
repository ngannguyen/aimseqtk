#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions to make Latex Tables and Text Tables
'''

import aimseqtk.lib.common as libcommon


def table_text_tab(f, objs, colfields):
    for (name, obj) in objs:
        f.write("%s" % name)
        for attr in colfields:
            if not obj or attr not in obj.getitems():
                f.write("\tNA")
            else:
                f.write("\t%.3f" % obj[attr])
                stdattr = "%s_std" % attr
                if stdattr in obj.getitems():
                    std = obj[stdattr]
                    if std:
                        f.write(" +/- %.3f" % std)
        f.write("\n")

def table_text_tab_list(f, objs, listattr):
    for (name, obj) in objs:
        if obj and listattr in obj.getitems():
            vals = ["%.3f" % v for v in obj[listattr]]
            f.write("%s\t%s\n" % (name, "\t".join(vals)))

def table_text(objs, outfile, colfields, listattr=None):
    # Rows = Samples, Columns = colfields
    f = open(outfile, 'w')
    f.write("#Sample\t%s\n" % "\t".join(colfields))
    if listattr:
        table_text_tab_list(f, objs, listattr)
    else:
        table_text_tab(f, objs, colfields)
    f.close()

def table_latex_tab(f, objs, colfields):
    for (name, obj) in objs:
        f.write("%s" % name.replace('_', '\_'))
        for attr in colfields:
            if not obj or attr not in obj.getitems():
                f.write(" & NA")
            else:
                f.write(" & %.3f" % obj[attr])
                stdattr = "%s_std" % attr
                if stdattr in obj.getitems():
                    std = obj[stdattr]
                    if std:
                        f.write(" $\pm$ %.3f" % std)
        f.write("\\\\\n")
        f.write("\\hline\n")

def table_latex_tab_list(f, objs, listattr):
    for (name, obj) in objs:
        if obj and listattr in obj.getitems():
            f.write("%s & " % name.replace('_', '\_'))
            vals = ["%.3f" % v for v in obj[listattr]]
            f.write(" & ".join(vals))
            f.write("\\\\\n")
            f.write("\\hline\n")

def table_latex(objs, outfile, colfields, caption='', label='',
                                                   listattr=None):
    f = open(outfile, 'w')
    libcommon.write_doc_start(f)
    colnames = ["Sample"] + colfields
    colnames = [c.replace('_', '\_') for c in colnames]
    libcommon.tab_header(f, colnames)
    if listattr:  # if want to print an obj attr which is a list
        table_latex_tab_list(f, objs, listattr)
    else:  # print all colfields attr
        table_latex_tab(f, objs, colfields)
    #caption = "Diversity colfields"
    libcommon.tab_closer(f, caption, label)
    libcommon.write_doc_end(f)
    f.close()

def table(name2obj, outfile, colfields, group2avr={}, group2names={},
          tex=False, keyattr='numclone', caption='', label='', islist=False):
    # rows=samples; cols=colfields;
    if not name2obj or not colfields:
        return
    # Sort the samples by group and/or size
    keyfunc = None
    listattr = None
    if islist or keyattr not in colfields:
        listattr = keyattr
    else:
        keyfunc = lambda item: item[0][keyattr]
    
    if group2avr:
        sortedobjs = libcommon.sort_objs_by_group(name2obj, group2names,
                                           True, group2avr, keyfunc=keyfunc)
    else:
        sortedobjs = libcommon.sort_dict_by_value(name2obj, keyfunc=keyfunc)

    if tex:
        table_latex(sortedobjs, outfile, colfields, caption, label, listattr)
    else:
        table_text(sortedobjs, outfile, colfields, listattr)

def matrix_table(names, rows, outfile):
    # Print matrix to tab-separated file
    assert len(names) == len(rows)
    f = open(outfile, 'w')
    f.write("Samples\t%s\n" % "\t".join(names))
    for i, name in enumerate(names):
        row = rows[i]
        assert len(names) == len(row)
        f.write("%s\t%s\n" % (name, "\t".join(["%.3f" % r for r in row])))
    f.close



