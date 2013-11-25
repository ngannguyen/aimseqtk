#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Draw rarefaction plots
'''

import os
import sys

import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties

import aimseqtk.lib.drawcommon as drawcommon


def draw_rarefaction(outfile, name2size2sampling, group2samples):
    # xaxis: sampling size; yaxis: index value +- std
    

