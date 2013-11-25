#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

'''
Common functions for plotting using Matplotlib
'''

import os
import sys


def get_colors_medium():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#377EB8", "#E31A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#1B9E77",
              "#FFFF33", "#A65628", "#CE1256"]

def get_colors_light():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#A6D7FE", "#FE8E8F", "#B8FEB5", "#F6BDFE", "#FEBF80", "#95FEDF",
              "#FFFFB3", "#D8885A", "#D7B5D8"]

def get_colors_dark():
    # blue, red, green, purple, orange, green-ish, yellow, brown, pink 
    colors = ["#275880", "#B30000", "#367A33", "#6A3772", "#B25900", "#136E53",
              "#B2B324", "#743C1C", "#900D3D"]

def get_markers():
    markers = ['o', 'd', '^', 'p', 'v', '*', 's', '+', 'x']
    return markers

def get_n_colors(n):
   import colorsys
   hsv_tuples = [(x*1.0/n, 0.5, 0.5) for x in xrange(n)]
   rgb_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples)
   return rgb_tuples

def getname2color(names):
    colors = get_colors_medium()
    if len(names) > len(colors):
        colors = get_n_colors(len(names))
    name2color = {}
    for i, color in enumerate(colors):
        name2color[names[i]] = color
    return name2color


