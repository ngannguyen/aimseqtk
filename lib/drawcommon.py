#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions for plotting using Matplotlib
'''

import os
import sys
from optparse import OptionGroup

import matplotlib
matplotlib.use('Agg')

def get_colors_medium():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#377EB8", "#E31A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77",
              "#FFFF33", "#A65628", "#CE1256"]
    return colors

def get_colors_light():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#A6D7FE", "#FE8E8F", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF",
              "#FFFFB3", "#D8885A", "#D7B5D8"]
    return colors

def get_colors_dark():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#275880", "#B30000", "#367A33", "#6A3772", "#B25900", "#136E53",
              "#B2B324", "#743C1C", "#900D3D"]
    return colors

def get_markers():
    markers = ['o', 'd', '^', 'p', 'v', '*', 's', '+', 'x']
    return markers

def get_n_colors(n):
   import colorsys
   hsv_tuples = [(x*1.0/n, 0.5, 0.5) for x in xrange(n)]
   rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
   return rgb_tuples

def get_name2color(names):
    colors = get_colors_medium()
    if len(names) > len(colors):
        colors = get_n_colors(len(names))
    name2color = {}
    for i, name in enumerate(names):
        name2color[name] = colors[i]
    return name2color

def get_name2color_wtgroup(names, name2group, group2names):
    name2color = {}
    if group2names:
        groups = group2names.keys()
        group2color = get_name2color(groups)
        for name, group in name2group.iteritems():
            if group in group2color:
                name2color[name] = group2color[group]
    else:
        name2color = get_name2color(names)
    return name2color

def get_name2marker(names):
    markers = get_markers()
    if len(names) > len(markers):
        markers = ['.'] * len(names)
    name2marker = {}
    for i, name in enumerate(names):
        name2marker[name] = markers[i]
    return name2marker

def get_name2marker_wtgroup(samplenames, group2names):
    markers = get_markers()
    name2marker = {}
    if group2names:
        for group, names in group2names.iteritems():
            if len(markers) < len(names):
                return {}
        for group, names in group2names.iteritems():
            for i, name in enumerate(names):
                name2marker[name] = markers[i]
    elif len(markers) >= len(samplenames):
        for i, name in enumerate(samplenames):
            name2marker[name] = markers[i]
    return name2marker

#======== Matplotlib related functions =========
def add_plot_options(parser):
    group = OptionGroup(parser, "Plot options")
    group.add_option('--makeplots', dest='makeplots', action='store_true',
                     default=False, help=('If specified, will produce plots.'
                     + ' Default=%default.'))
    group.add_option('--plotformat', dest='plotformat', default='pdf',
                     help=('Plot output format [pdf|png|eps|all]. ' +
                           'Default=%default'))
    group.add_option('--dpi', dest='dpi', default=300, type='int',
                     help='Dots per inch. Default=%default')
    parser.add_option_group(group)

def check_plot_options(parser, options):
    if options.dpi < 72:
        parser.error('dpi must >= than screen res, 72. Got: %d' % options.dpi)
    if options.plotformat not in ('pdf', 'png', 'eps', 'all'):
        parser.error('Unrecognized plot format: %s.' % options.plotformat +
                     ' Choose one from: pdf png eps all.')

def init_image(width, height, outformat, outname, dpi):
    """
    init_image takes a width and height and returns
    both a fig and pdf object.
    """
    #matplotlib.use('PDF')
    import matplotlib.backends.backend_pdf as pltback
    import matplotlib.pyplot as plt
    pdf = None
    if outformat == 'pdf' or outformat == 'all':
        pdf = pltback.PdfPages(outname + '.pdf')
    fig = plt.figure(figsize=(width, height), dpi=dpi, facecolor='w')
    return (fig, pdf)

def write_image(fig, pdf, outformat, outname, dpi):
    if outformat == 'pdf':
        fig.savefig(pdf, format='pdf')
        pdf.close()
    elif outformat == 'png':
        fig.savefig(outname + '.png', format='png', dpi=dpi)
    elif outformat == 'eps':
        fig.savefig(outname + '.eps', format='eps')
    elif outformat == 'all':
        fig.savefig(pdf, format='pdf')
        pdf.close()
        fig.savefig(outname + '.png', format='png', dpi=dpi)
        fig.savefig(outname + '.eps', format='eps')

def set_axes(fig):                                                                     
    # Set single axes
    return fig.add_axes([0.12, 0.15, 0.85, 0.75])                                      

def get_axes(outfile="plot", w=10.0, h=8.0, outfmt='pdf', dpi=300):
    fig, pdf = init_image(w, h, outfmt, outfile, dpi)
    axes = fig.add_axes([0.12, 0.15, 0.85, 0.75])                                      
    return axes, fig, pdf

def set_axes2(fig, numsam, samples_per_plot):
    # Set multiple axes, depending on the number of samples_per_plot
    # and total Samples
    numaxes = numsam / samples_per_plot 
    if numsam % samples_per_plot != 0:
        numaxes += 1
    
    axes_list = []
    axleft = 0.12
    axright = 0.97
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.93
    axheight = axtop - axbottom
    margin = 0.1

    h = (axheight - (margin * (numaxes - 1)))/float(numaxes)

    bottom = axtop - h
    for i in xrange(numaxes):#add axes from top down
        axes_list.append(fig.add_axes([axleft, bottom, axwidth, h]))
        bottom = bottom - (margin + h)
    return axes_list

def set_axes3(fig, range1, range2):
    axleft = 0.12
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.05

    h1 = (axheight - margin)*(range1/(range1 + range2))
    h2 = axheight - margin - h1

    ax2 = fig.add_axes([axleft, axbottom, axwidth, h2])
    ax = fig.add_axes([axleft, axbottom + h2 + margin, axwidth, h1])
    return ax, ax2

def draw_discontinue_sign(top_axes, bottom_axes, top, bottom):
    d = 2
    top_axes.plot((-0.6, -0.4), (top + d, top - d), color="k", clip_on=False)
    bottom_axes.plot((-0.6,-0.4), (bottom +d, bottom -d), color="k",
                     clip_on=False)

def edit_spine2(top_axes, bottom_axes):
    # Set spines for plot with discontinous sign
    top_axes.spines['bottom'].set_visible(False)
    top_axes.spines['top'].set_visible(False)
    top_axes.spines['right'].set_visible(False)
    top_axes.yaxis.set_ticks_position('left')
    top_axes.xaxis.set_ticks_position('none')

    bottom_axes.spines['top'].set_visible(False)
    bottom_axes.spines['right'].set_visible(False)
    bottom_axes.xaxis.tick_bottom()
    bottom_axes.yaxis.set_ticks_position('left')

def edit_spine(axes):                                                                 
    for loc, spine in axes.spines.iteritems():
        if loc in ['left', 'bottom']:                                                 
            spine.set_position(('outward', 10))                                       
        elif loc in ['right', 'top']:
            spine.set_color('none')
        else:
            raise ValueError('Unknown spine location %s\n' % loc)                     

def set_grid(axes):
    axes.xaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)
    axes.yaxis.grid(b=True, color='#3F3F3F', linestyle='-', linewidth=0.05)

def set_xticks(axes, data, labels):
    axes.xaxis.set_ticks(data)
    axes.xaxis.set_ticklabels(labels)

def set_yticks(axes, data, labels):
    axes.yaxis.set_ticks(data)
    axes.yaxis.set_ticklabels(labels)

def set_labels(axes, title=None, xlabel=None, ylabel=None):
    if title:
        axes.set_title(title, size='xx-large', weight='bold')
    if xlabel:
        axes.set_xlabel(xlabel, size='x-large', weight='bold')
    if ylabel:
        axes.set_ylabel(ylabel, size='x-large', weight='bold')

def set_ticks(axes):
    axes.xaxis.set_ticks_position('bottom')                                           
    axes.yaxis.set_ticks_position('left')                                             
    minorLocator = LogLocator(base=10, subs = range(1, 10))                           
    axes.xaxis.set_minor_locator(minorLocator)                                        

def adjust_xticklabels(axes, size='small', weight='bold', rotation=None):
    for label in axes.get_xticklabels():
        label.set_fontsize(size)
        label.set_fontweight(weight)
        if rotation:
            label.set_rotation(rotation)

def adjust_yticklabels(axes, size='small', weight='bold', rotation=None):
    for label in axes.get_yticklabels():
        label.set_fontsize(size)
        label.set_fontweight(weight)
        if rotation:
            label.set_rotation(rotation)

def adjust_ticklabels(axes, size='small', weight='bold', xrotation=None,
                                                         yrotation=None):
    adjust_xticklabels(axes, size, weight, xrotation)
    adjust_yticklabels(axes, size, weight, yrotation)

def set_legend(axes, lines, linenames):
    legend = axes.legend(lines, linenames, numpoints=1, loc='best', ncol=1)
    legend._drawFrame = False

def bihist(y1, y2, axes, bins, orientation, color=None):
    # Top hist
    n1, bins1, patch1 = axes.hist(y1, bins=bins, orientation=orientation, 
                                  color=color)
    # Bottom hist
    n2, bins2, patch2 = axes.hist(y2, bins=bins, orientation=orientation,
                                  color=color)
    #set ymax:
    if orientation == 'vertical':
        ymax = max([i.get_height() for i in patch1])
        for i in patch2:
            i.set_height(-i.get_height())
        ymin = min([i.get_height() for i in patch2])
    elif orientation == 'horizontal':
        ymax = max([i.get_width() for i in patch1])
        for i in patch2:
            i.set_width(-i.get_width())
        ymin = min([i.get_width() for i in patch2])
    #axes.set_ylim(ymin*1.1, ymax*1.1)
    return ymin, ymax

def draw_heatmap(rownames, colnames, rows, outfile):
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    gplots = importr('gplots')
    grdevices = importr('grDevices')
    # convert rows into r matrix
    mvec = []
    for row in rows:
        mvec.extend(row)
    rvec = ro.FloatVector(mvec)
    rrownames = ro.StrVector(rownames)
    rcolnames = ro.StrVector(colnames)
    rmatrix = ro.r['matrix'](rvec, nrow=len(rows), byrow=True)

    print "outfile %s" % outfile
    grdevices.pdf(file=outfile)
    try:
        #gplots.heatmap_2(rmatrix)
        gplots.heatmap_2(rmatrix, Rowv=True, Colv=True, labRow=rrownames,
                                                        labCol=rcolnames)
    except:
        print rmatrix
        print rrownames
        print rcolnames
        #raise RuntimeError("Heatmap errors")
    grdevices.dev_off()


