#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 08:18:40 2020

this is the calculator for delta
read Das & Pappu for details on the definition of each term and the meaning
https://doi.org/10.1073/pnas.1304749110

@author: yiwenchen
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import seaborn as sns


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 20})
sns.set(style='white')

def make_cmap(colors, position=None, bit=False, norm=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # Draw the heatmap with the mask and correct aspect ratio
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    if norm != False:
        norm = mpl.colors.BoundaryNorm(np.linspace(0,1,len(colors)+1), cmap.N)
    return cmap, norm


def sigma(seq):
    f_pos = seq.count('K')*1.0/len(seq)
    f_neg = seq.count('E')*1.0/len(seq)
    out = (f_pos-f_neg)**2/(f_pos+f_neg)
    return out

def delta(sigma_array):
    avg = sigma_array.mean()
    out = 0
    for i in sigma_array:
        out += (i-avg)**2
    out = out*1.0/len(sigma_array)
    return out
    
def seq_to_num(seqs):
    matrix = np.zeros((len(seqs), 50))
    for i, seq_i in enumerate(seqs):
        for j, aa in enumerate(seq_i):
            matrix[i,j] = 1 if aa == 'K' else 0.5
    return matrix


# this function is used to generate disorder score in the form of a compressed heatmap 
# it is adpoted from the two_direction_heatmap in CorrelomicsHeatmap_Final.py script
def par_strip(matrix, out_path, out_fname, mask='None', \
                          basecolor = None, ylabels = True, map_square = True, \
                      figsize = (12,3), yaxis_rot = 0, ax = None, cbar= None):
    # here the blue is for "-" charged, red for "+" charged
    cmap, norm = make_cmap([(0.95703125,0.95703125,0.957031255),(1/256.,170/256.,255/256.),\
                            (255/256.,1/256.,81/256.)], norm=True)
    if mask == 'None': 
        mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    elif mask == True:
        mask = matrix.isna()
    # specify kws for cbar
    cbar_kws_dict = {"shrink":0.5, "location":'top', "use_gridspec":False} if cbar == 'horizontal' else {"shrink":0.5}
    if ax == None:
        f, ax = plt.subplots(figsize=figsize)
        save = True
    else:
        save = False
    sns.set(font_scale = 2)
    sns.set(style='white')
    if norm != False:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
                square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
    else:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin = 0, vmax=1,\
                square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
    ax.set_facecolor('lightgrey')
    ax.hlines(np.arange(1, matrix.shape[0]), *ax.get_xlim(), color= 'white')
    matrix_map.tick_params(labelsize=12)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    if save == True:
        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    else:             
        ax.set_ylabel('')
        ax.set_aspect(matrix.shape[1]/200.*10)
    plt.close()
    return ax



sv1 = 'EK'*25
sv2 = 'EEEKKK'*8+'EK'
sv3 = 'KEKKKEKKEEKKEEKEKEKEKEEKKKEEKEKEKEKKKEEKEKEEKKEEEE'
sv6 = 'EEEKKEKKEEKEEKKEKKEKEEEKKKEKEEKKEEEKKKEKEEEEKKKKEK'
sv8 = 'KKKKEEEE'*6+'KE'
sv12 = 'EKKEEEEEEKEKKEEEEKEKEKKEKEEKEKKEKKKEKKEEEKEKKKKEKK'
sv15 = 'KKEKKEKKKEKKEKKEEEKEKEKKEKKKKEKEKKEEEEEEEEKEEKKEEE'
sv19 = 'EEEEEKKKKK'*5
sv28 = 'E'+'K'*21+'E'*18+'KKEEEEEKEK'
sv29 = 'KEEEEK'+'E'*21+'K'*23
sv30 = 'E'*25+'K'*25


########## Section 1: get a list of sequences, then compute the delta
seq_list = [sv1, sv2, sv3, sv6, sv8, sv12,sv30, sv15, sv19, sv29, sv28]
delta_list = np.zeros(len(seq_list))
for i, sv_i in enumerate(seq_list):
    delta_i = 0
    for j, g_j in enumerate([5,6]):
        N_blob = len(sv_i)-g_j
        sigma_array = np.zeros(N_blob)
        for z in range(N_blob):
            seq_z = sv_i[z:z+g_j]
            sigma_array[z] = sigma(seq_z)
        delta_i += delta(sigma_array)
#        print delta(sigma_array)
    delta_list[i] = delta_i/(j+1)
# part 2. generate example plot to illustrate delta difference among peptide with different charge spacings
seqs = [sv1, sv2, sv3, sv6, sv8, sv12,sv30, sv15, sv19, sv29, sv28]
deltas = delta_list
# deltas =[2.40741243e-35, 2.40741243e-35, 7.37327946e-03, 1.26542130e-02,
#        1.39510864e-02, 4.64071576e-02, 5.94835537e-02, 7.28683806e-02,
#        7.85428466e-02, 9.27587928e-02, 1.34287943e-01]
matrix = seq_to_num(seqs)
df_delta = pd.DataFrame(data = matrix, index = np.round(deltas,4))
fpath_out = 'Properties/delta'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
fname_out = 'delta_strip.pdf'
ax = par_strip(df_delta, fpath_out, fname_out)   
