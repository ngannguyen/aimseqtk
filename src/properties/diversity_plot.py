#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw diversity plots
xaxis = groups
yaxis = diversity index values
'''

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties

import aimseqtk.lib.drawcommon as drawcommon


def draw_diversity_plot(group2names, name2obj, attr, outfile, outfmt):
     



