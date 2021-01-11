#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 20:38:02 2019

@author: yiwenchen
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from skimage import io
import pylab
import os
import pandas as pd

#def rgb(minimum, maximum, value):
#    minimum, maximum = float(minimum), float(maximum)
#    ratio = 2 * (value-minimum) / (maximum - minimum)
#    b = int(max(0, 255*(1 - ratio)))
#    r = int(max(0, 255*(ratio - 1)))
#    g = 255 - b - r
#    return [r, g, b]
#
##rgb_min = rgb(-0.5, 0.5, -0.5)
##rgb_middle = rgb(-0.5, 0.5, 0)
##rgb_max = rgb(-0.5, 0.5, 0.5)
##
##palette = np.array([rgb_min, rgb_middle, rgb_max], dtype = np.uint8)
##indices = np.array([[0, 1, 2]])
##io.imshow(palette[indices])


def linear_rescale(x, vmin, vmax, original_scale = 1):
    if x > vmax:
        x_rescale = 1
    elif x < vmin:
        x_rescale = 0
    else:
        x_rescale = (x - vmin)*1. /(vmax-vmin) * original_scale
    return x_rescale


def custom_color_hex_old(x, vmin, vmax):
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"]) 
    return mpl.colors.to_hex(cmap(linear_rescale(x, vmin, vmax)))

def custom_color_hex(x, vmin, vmax):
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"]) 
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    return mpl.colors.to_hex(cmap(norm(x)))


def sub_table(df, tissues, var):
    df_new = df.loc[0:,[i for i in df.columns if var in i]]
    df_new.columns = df_new.iloc[0]
    df_new = df_new[tissues]
    return df_new

# this create a simple standalone color bar
def color_bar(vmin, vmax, cmap=None):
    # Make a figure and axes with dimensions as desired.
    fig, ax = plt.subplots(figsize=(8, 2))  
    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    if cmap == None: cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"]) 
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
    cb.set_label('example colorbar')
    fig.savefig('colorbar_%s_%s.pdf'%(vmin, vmax))


### useful reference link for the section below:
### https://matplotlib.org/tutorials/colors/colorbar_only.html
### https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
## Creat a colormap with a gradient from blue to white between -0.5 and 0,\
### then white to yellow from 0 to 0.5.
#fig, ax = plt.subplots(figsize=(2, 12))
#fig.subplots_adjust(bottom=0.5)
#cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"]) #this is set to blue when x = 0, white when x =0.5, yellow when x = 1
## if I want the actual colorbar as colors = [(-0.5, "blue"), (0, "white"), (0.5, "yellow")], the run the next 3 lines
#norm = mpl.colors.Normalize(vmin = -2, vmax = 2)
#cb1 = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')
#fig.savefig('colobarv.pdf')
#fig.show()

### so if I want to know the hex value of logFC = 0.5, then type mpl.colors.to_hex(cmap(1)),\
### note that I used 1 because logFC 0.5 is scaled to 1 in the colorbar, alternatively do this\
### mpl.colors.to_hex(cmap(linear_rescale(x, -0.5, 0.5))) where x = 0.5
#
# read the logFC value of a table, then output all the hex value based on a specified colorscale
# create a linear colormap of interest
#cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])    
#cmap_vmin = -0.5; cmap_vmax = 0.5
#f_table = 'YC_BGLPMHMST_AllChaperones_OvY.csv'
#tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin']#,'Testis']
#table_path = 'GeneSets/MS/Chaperones'
#var = 'rankstats'
#df_table = pd.read_table(os.path.join(table_path, f_table),sep=',')
#df_table_filter = sub_table(df_table, tissues, var)
##df_table = pd.read_table(os.path.join(table_path, f_table), sep = '\t')
#for i in tissues:
#    df_table_filter ['%s_hex'%i] = df_table_filter[i].apply(lambda x: mpl.colors.to_hex(cmap(linear_rescale(x, cmap_vmin, cmap_vmax))))
#df_table_filter.to_csv(os.path.join(table_path, f_table.replace('.txt', '_%s_hex.csv'%var)), index = False)



###### Updated on 2020.01.11 #######
cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])    
cmap_vmin = -2; cmap_vmax = 2
#### part 1. get the hex color on proteasome
f_table = 'Proteasome_All_OvY_AGG_flat.csv'
tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis']
table_path = 'GeneSets/MS/Proteasome_AGG'
var = 'rankstats'
df_table = pd.read_csv(os.path.join(table_path, f_table))
df_table['%s_hex'%var] = df_table[var].apply(lambda x: mpl.colors.to_hex(cmap(linear_rescale(x, cmap_vmin, cmap_vmax))))
df_table.to_csv(os.path.join(table_path, f_table.replace('.csv', '_%s_hex.csv'%var)), index = False)
#### part 2. get the hex color on TRiC complex
f_table = 'Chaperonines_All_OvY_AGG_flat.csv'
tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis']
table_path = 'GeneSets/MS/Chaperones_AGG'
var = 'rankstats'
df_table = pd.read_csv(os.path.join(table_path, f_table))
df_table['%s_hex'%var] = df_table[var].apply(lambda x: mpl.colors.to_hex(cmap(linear_rescale(x, cmap_vmin, cmap_vmax))))
df_table.to_csv(os.path.join(table_path, f_table.replace('.csv', '_%s_hex.csv'%var)), index = False)