#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 07:41:02 2019

this Diso__Heatmap is used to visualize the disorder results using color scale and stack multiple sequences at the same time

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

# function used to parse a parameter by delimiter for downstream applications (i.e. parse the Diso Score by "-"
# to get the individual score for each amino acid residue).
def coln_to_matrix(df, identifier, par_name, par_delimiter, in_list, order = None):
    df_plot = df[[identifier,par_name]][df[identifier].isin(in_list)]
    df_split = df_plot[par_name].str.split(pat=par_delimiter,expand = True)
    df_split.set_index(df_plot.loc[df_split.index.tolist(),identifier], inplace=True)
    if order == None:
        order = df_split.index.tolist()
    df_out = df_split.loc[order,:].astype('float64')
    df_out.dropna(axis = 'columns', how = 'all', inplace = True)
    df_mask = df_out.isna()
    return df_out, df_mask


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


# function used to make customized cmap for different needs
def pick_cmap(mtype, basecolor=None):
    if mtype == 'Diso Score':
        if basecolor == None:
#            cmap, norm = make_cmap([(1,1,1),(0,0,1)], position=[0,1], bit=False) # this is going from white to blue
            cmap, norm = make_cmap([(1,1,1),(0.9296875, 0.5078125, 0.9296875)], position=[0,1], bit=False) # this is going from white to purple
        else:
            cmap, norm = make_cmap([(1,1,1),mpl.colors.to_rgb(basecolor)],  position=[0,1], bit=False)
    elif mtype == 'local uversky':
        cmap, norm = make_cmap([(1,1,1),mpl.colors.to_rgb('grey')],  position=[0,1], bit=False)
    elif mtype == 'hybrid':
        cmap, norm = make_cmap([(0.95703125,0.95703125,0.957031255), (1/256.,170/256.,255/256.),\
                                (255/256.,1/256.,81/256.),(255/256.,213/256.,1/256.)], norm = True)
#        cmap, norm = make_cmap([(0.95703125,0.95703125,0.957031255), (0.984375,0.55078125,0.3828125),\
#                                (0.55078125,0.625,0.79296875), (0.3984375,0.7578125,0.64453125)], norm = True)
    elif mtype == 'NCPR':
        cmap, norm = make_cmap([(0.95703125,0.95703125,0.957031255),(1/256.,170/256.,255/256.),\
                                (255/256.,1/256.,81/256.)], norm=True)
    elif mtype == 'hydropathy':
        cmap, norm = make_cmap([(0.95703125,0.95703125,0.957031255),(255/256.,213/256.,1/256.)], norm = True)
    else:
        aa_value_dict, cmap, norm = aa_color(mtype)
    return cmap, norm


# this file is used to create a custom colormap where different AA are shown in different colors
# note that the decimal points associated with each type is the basis of the order of position
def aa_color(aa_type, type_color_dict=None):
    if type_color_dict == None:
        if aa_type == 'charge':
            type_color_dict = {'non-polar aliphatic':['lightgrey',0], 'non-polar aromatic':['lightgrey',0.2],\
                       'polar negative':['blue',0.4], 'polar positive':['red',0.6], 'polar neutral':['lightgrey',0.8],\
                       'other':['lightgrey',1]}
        elif aa_type == 'aromatic':
            type_color_dict = {'non-polar aliphatic':['lightgrey',0], 'non-polar aromatic':['green',0.2],\
                       'polar negative':['lightgrey',0.4], 'polar positive':['lightgrey',0.6], 'polar neutral':['lightgrey',0.8],\
                       'other':['lightgrey',1]}
        else:
            type_color_dict = {'non-polar aliphatic':['lightgrey',0], 'non-polar aromatic':['green',0.2],\
                       'polar negative':['blue',0.4], 'polar positive':['red',0.6], 'polar neutral':['lightgrey',0.8],\
                       'other':['lightgrey',1]}                  
    aa_type_dict = {'non-polar aliphatic':['A','V','L','I','M'], 'non-polar aromatic':['F','Y','W'],\
                       'polar negative':['D','E'], 'polar positive':['H','K','R'], 'polar neutral':['S','T','N','Q'],\
                       'other':['P','G','C','B','Z']}
    aa_value_dict = dict()
    df_cp = pd.DataFrame.from_dict(type_color_dict, orient='index', columns = ['Color','Position'])
    for type_i in type_color_dict.keys():
        color_i, position_i = type_color_dict[type_i]
        for aa_i in aa_type_dict[type_i]:
            aa_value_dict[aa_i] = position_i
    df_cp.sort_values(by=['Position'], ascending = True, inplace=True)
    aa_cmap, aa_norm = make_cmap([mpl.colors.to_rgb(c) for c in df_cp.Color.tolist()], df_cp.Position.tolist(), norm = True)
    return aa_value_dict, aa_cmap, aa_norm
    

# this function is used to convert a protein sequence to decimal chain to visualize the AA residues
def aa_to_numeric(aa, aa_value_dict):
    aa_num = aa
    for aa_i in aa_value_dict.keys():
        aa_num = aa_num.replace(aa_i, '-%s'%aa_value_dict[aa_i])
    return aa_num[1:]


# this function is used to generate disorder score in the form of a compressed heatmap 
# it is adpoted from the two_direction_heatmap in CorrelomicsHeatmap_Final.py script
def par_strip(matrix, out_path, out_fname, mtype, cutoff = 0.05, mask='None', \
                          basecolor = None, ylabels = True, map_square = True, \
                      figsize = (12,3), yaxis_rot = 0, ax = None, cbar= None):
    cmap, norm = pick_cmap(mtype=mtype, basecolor=basecolor)
    matrix.drop_duplicates(inplace=True)
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
    matrix_map.tick_params(labelsize=12)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    if save == True:
        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    else:             
        ax.set_ylabel('')
        ax.set_aspect(matrix.shape[1]/200.*10)
    plt.close()
    return ax


## this function is used to generate disorder score in the form of a compressed heatmap 
## it is adpoted from the two_direction_heatmap in CorrelomicsHeatmap_Final.py script
#def par_strip(matrix, out_path, out_fname, mtype, cutoff = 0.05, mask='None', \
#                          basecolor = None, ylabels = True, map_square = True, \
#                      figsize = (12,3), yaxis_rot = 0, ax = None, cbar= None):
#    cmap, norm = pick_cmap(mtype=mtype, basecolor=basecolor)
#    matrix.drop_duplicates(inplace=True)
#    if mask == 'None': 
#        mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
#    elif mask == True:
#        mask = matrix.isna()
#    # specify kws for cbar
#    if cbar == 'horizontal':
#        cbar_kws_dict = {"shrink":0.5, "location":'top', "use_gridspec":False}
#    else:
#        cbar_kws_dict = {"shrink":0.5}
#    if ax == None:
#        f, ax = plt.subplots(figsize=figsize)
#        sns.set(font_scale = 2)
#        if norm != False:
#            matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
#                    square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
#        else:
#            matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin = 0, vmax=1,\
#                    square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
#        ax.set_facecolor('lightgrey')
#        # without linewidth=.5, the figure looks much better
#        #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
#        #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
#        matrix_map.tick_params(labelsize=12)
#        matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
#        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
#        plt.close()
#    else:
#        sns.set(font_scale = 2)
#        sns.set(style='white')
#        if norm != False:
#            matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
#                    square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
#        else:
#            matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin = 0, vmax=1,\
#                    square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels) 
##            matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin = 0, vmax=1,\
##                    square=map_square, ax = ax, cbar_kws=cbar_kws_dict, yticklabels = ylabels, xticklabels=False)              
#        ax.set_ylabel('')
#        ax.set_aspect(matrix.shape[1]/200.*10)
#        ax.set_facecolor('lightgrey')
#        # without linewidth=.5, the figure looks much better
#        matrix_map.tick_params(labelsize=12)
#        matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
#    return ax

        
# plot disorder trace
def diso_trace(df, out_path, out_fname, identifier, key, par='Diso Score', par_delimiter = '-', ax = False):
    df_trace = df[[identifier,par]].drop_duplicates()
    diso_str = df_trace[par][df_trace[identifier]==key].values
    diso_trace_array = np.asarray(diso_str[0].split(par_delimiter)).astype(float)
    x = range(len(diso_trace_array))
    sns.set(style="white")
    if ax == False:
        f, ax = plt.subplots(figsize = (12,2))
        new_fig = True
    else:
        new_fig = False
    if 'NCPR' in par:
        ax.fill_between(x, 0, diso_trace_array, where=0 >=diso_trace_array, facecolor='blue', edgecolor='blue')#,alpha= 0.5)
        ax.fill_between(x, 0, diso_trace_array, where=0 <=diso_trace_array, facecolor='red', edgecolor='red')#,alpha= 0.5)
        ax.set_ylim(-1,1)            
    elif 'hydropathy' in par:
        ax.plot(x, diso_trace_array,'lightgrey')
        ax.plot(x, [0.458484]*len(x), 'k--')
        ax.fill_between(x, 0, diso_trace_array, where=0 <=diso_trace_array, facecolor='lightgrey')
        ax.fill_between(x, 0.458484, diso_trace_array, where=0.458484 <=diso_trace_array, facecolor='red', edgecolor='red')
        ax.set_ylim(0,1)
    else:
        ax.fill_between(x, 0, diso_trace_array, where=0 <=diso_trace_array, facecolor='lightgrey')
        ax.set_ylim(0,1)
        ax.plot(x, diso_trace_array,'lightgrey')
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(6))
    ax.set_ylabel(par)
    if new_fig == True:
        f.savefig(os.path.join(out_path, out_fname))
        plt.close()
    else:
        return ax


# local uverskyplot factor
def uversky_local(NCPR, hydropathy):
    #mean_hydropathy > (MeanNetCharge+1.151)/2.785 ==>folded # check reference at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373528/
    cutoff = (abs(NCPR)+1.151)/2.785
    folded = np.greater_equal(hydropathy, cutoff)
    return folded


# get a local_uversky plot
def uversky_local_trace(df, out_path, out_fname, identifier, key, par_delimiter = ',', ax=False):
    NCPR_coln = 'NCPR_blobLen5'
    hydropathy_coln = 'hydropathy_blobLen5'
    df_key = df[[identifier, NCPR_coln, hydropathy_coln]].drop_duplicates()
    NCPR_array = np.asarray(df_key[NCPR_coln][df_key[identifier]==key].values[0].split(par_delimiter)).astype(float)
    hydropathy_array = np.asarray(df_key[hydropathy_coln][df_key[identifier]==key].values[0].split(par_delimiter)).astype(float)
    x = range(len(NCPR_array))
    folded = uversky_local(NCPR_array, hydropathy_array)
    sns.set(style='white')
    if ax == False:
        f, ax = plt.subplots(figsize = (12,2))
        new_fig = True
    else:
        new_fig = False
#    ax.set_ylim(0,1)
#    ax.plot(x, folded,'lightgrey')
#    ax.fill_between(x, 0, folded, where=0 <=folded, facecolor='lightgrey')
#    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(6))
#    ax.set_ylabel('folded (uversky)')
    df_folded = pd.DataFrame(np.array([x, ~folded]), index = ['index', 'unfolded'], columns = range(len(folded)))
    ax = par_strip(df_folded.iloc[1:], out_path, out_fname, 'local uversky', cutoff = 0.05, mask='None', ax = ax, cbar='horizontal')
    if new_fig == True:
        f.savefig(os.path.join(out_path, out_fname))
        plt.close()
    else:
        return ax
    

def NCPR_hydropathy_hybrid(df, out_path, out_fname, identifier, key, ptype='hybrid', par_delimiter = ',', ax=False):
    NCPR_coln = 'NCPR_blobLen5'
    hydropathy_coln = 'hydropathy_blobLen5'
    df_key = df[[identifier, NCPR_coln, hydropathy_coln]].drop_duplicates()
    NCPR_array = np.asarray(df_key[NCPR_coln][df_key[identifier]==key].values[0].split(par_delimiter)).astype(float)
    hydropathy_array = np.asarray(df_key[hydropathy_coln][df_key[identifier]==key].values[0].split(par_delimiter)).astype(float)
    df_out = pd.DataFrame({'NCPR':NCPR_array, 'hydropathy':hydropathy_array})
    df_out['plot'] = 0 # for NCPR is positive (0.5), negative (0)
    if ptype == 'hydropathy':
        df_out.loc[df_out['hydropathy']>=0.5,'plot']= 1  
    else:
        df_out.loc[df_out['NCPR']>0,'plot']= 0.7 ############### "+" charge
        df_out.loc[df_out['NCPR']<0,'plot']= 0.4 ############### "-" charge
        if ptype == 'hybrid':
            df_out.loc[df_out['hydropathy']>=0.5,'plot']= 1 # hydropathy above 0.5 will be plotted
    sns.set(style='white')
    if ax == False:
        f, ax = plt.subplots(figsize = (12,2))
        new_fig = True
    else:
        new_fig = False
#    ax.set_ylim(0,1)
#    ax.plot(x, folded,'lightgrey')
#    ax.fill_between(x, 0, folded, where=0 <=folded, facecolor='lightgrey')
#    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(6))
#    ax.set_ylabel('folded (uversky)'
    df_strip = pd.DataFrame(np.array([df_out.index, df_out['plot']]), index = ['index', ptype], columns = df_out.index)
    ax = par_strip(df_strip.iloc[1:], out_path, out_fname, ptype, ax = ax, cbar='horizontal')
    if new_fig == True:
        f.savefig(os.path.join(out_path, out_fname))
        plt.close()
    else:
        return ax         
    

# function use to generate a complied set of properties for a given protein
def par_compile(df, out_path, out_fname, identifier, key, pars = ['Diso Score', 'charge']):
    aa_length = df['AA Length'][df[identifier]==key].iloc[0]
    f, ax = plt.subplots(nrows=len(pars), ncols=1, figsize=(0.05*aa_length,3*len(pars)), sharex = True)
    for i, par_i in enumerate(pars):
        if (par_i == 'diso trace'):
            ax[i] = diso_trace(df, out_path, out_fname, identifier, key, par = 'Diso Score', ax = ax[i])
        elif ('blobLen' in par_i):
            ax[i] = diso_trace(df, out_path, out_fname, identifier, key, par = par_i, par_delimiter = ',', ax = ax[i])
        elif (par_i == 'charge') or (par_i == 'aromatic'):
            aa_value_dict, aa_cmap, aa_cmap_norm = aa_color(par_i)
            df_plot = df[[identifier,'Seq']][df[identifier]==key]
            df_plot['Seq Numeric'] = df_plot['Seq'].apply(lambda x: aa_to_numeric(x, aa_value_dict))
            par_delimiter='-'
            df_out, df_mask = coln_to_matrix(df_plot, identifier, 'Seq Numeric', par_delimiter, [key], order = [key])
            ax[i] = par_strip(df_out, out_path, out_fname, par_i, ax=ax[i], cbar = 'horizontal')      
        elif (par_i == 'local uversky'):
            ax[i] = uversky_local_trace(df, out_path, out_fname, identifier, key, par_delimiter = ',', ax=ax[i])
        elif (par_i == 'NCPR hydropathy hybrid'):
            ax[i] = NCPR_hydropathy_hybrid(df, out_path, out_fname, identifier, key, ptype='hybrid', ax=ax[i])
        elif (par_i == 'NCPR strip'):
            ax[i] = NCPR_hydropathy_hybrid(df, out_path, out_fname, identifier, key, ptype='NCPR', ax=ax[i])
        elif (par_i == 'hydropathy strip'):
            ax[i] = NCPR_hydropathy_hybrid(df, out_path, out_fname, identifier, key, ptype='hydropathy', ax=ax[i])
        else:
            if par_i == 'Diso Score': par_delimiter='-'
            df_out, df_mask = coln_to_matrix(df, identifier, par_i, par_delimiter, [key], order = [key])
            ax[i] = par_strip(df_out, out_path, out_fname, par_i, ax=ax[i], cbar = 'horizontal')
    f.savefig(os.path.join(out_path, out_fname))
    plt.close(f)

    
###### Output properties of proteins that are either AGG or Prop hits
fpath_out = 'Properties/AA_example'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
###### part 0. run this section to initialize the input folder.
if os.path.exists(os.path.join('Properties','AllDetected_metrics.csv')) == False:
    df_metrics = pd.read_csv(os.path.join('Properties','Killifish_DISOLLRCider.csv'))
    df_hit = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
    df_hit_nodup = df_hit[['N. furzeri Protein Id','N. furzeri Final Symbol']].drop_duplicates()
    df_MS = df_hit_nodup.merge(df_metrics, on = 'N. furzeri Protein Id', how='left')
    df_MS.to_csv(os.path.join('Properties','AllDetected_metrics.csv'),index=False)
else:
    df_MS = pd.read_csv(os.path.join('Properties','AllDetected_metrics.csv'))
###### part 1. run this to visualize sequence features of individual proteins
df_MS = pd.read_csv(os.path.join('Properties','AllDetected_metrics.csv'))
#pid = 'XP_015821815.1'; id_key = 'N. furzeri Protein Id' #casq2 (1 of 2), muscle, skin AGG hit
#pid = 'XP_015822974.1'; id_key = 'N. furzeri Protein Id' #bsn, brain agg hit
#pid = 'XP_015829630.1'; id_key = 'N. furzeri Protein Id' #sfpq, brain agg hit
pid = 'XP_015809881.1'; id_key  = 'N. furzeri Protein Id' #ddx5, brain prop hit
#pid = 'XP_015804619.1'; id_key  = 'N. furzeri Protein Id' #trio, brain AGG and Prop hit
#pid = 'XP_015827700.1'; id_key  = 'N. furzeri Protein Id' #nono, brain AGG and Prop hit
#pid = 'XP_015832600.1'; id_key  = 'N. furzeri Protein Id' #cct4, liver and skin Prop hit high aliphatic frac
#pid = 'XP_015811604.1'; id_key = 'N. furzeri Protein Id' #pex19, liver Prop hit
#pid = 'XP_015827480.1'; id_key = 'N. furzeri Protein Id' #ranbp1, gut AGG hit
#pid = 'XP_015799355.1 '; id_key = 'N. furzeri Protein Id' #hnrnpd, testis Prop hit
#pid = 'XP_015825944.1'; id_key = 'N. furzeri Protein Id' #SDSL testis Prop hit
#pid = 'XP_015801955.1'; id_key = 'N. furzeri Protein Id' #nop58
#pid = 'XP_015819799.1'; id_key = 'N. furzeri Protein Id' #slc25a22
out_name = '%s_biophysical.pdf'%df_MS['N. furzeri Final Symbol'][df_MS[id_key]==pid].unique()[0]
par_compile(df_MS, fpath_out, out_name, id_key, pid,\
            pars=['charge','aromatic','Diso Score', 'diso trace', 'NCPR_blobLen5',\
                  'hydropathy_blobLen5', 'FCR_blobLen5','sigma_blobLen5','NCPR hydropathy hybrid', 'NCPR strip', 'hydropathy strip'])
