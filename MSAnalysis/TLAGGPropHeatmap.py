#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:41:54 2019

this is used to generate heatmap for the correlomics/enrichment differences among different parameters
this is adapted from the CorrelomicsHeatmap_Final.py


@author: yiwenchen
"""

from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import itertools


mpl.rcParams['pdf.fonttype'] = 42

sns.set(style="white")


def two_scale_heatmap(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    f, ax = plt.subplots(figsize=figsize)
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(0.99609375, 0.625, 0.4765625), (0.859375, 0.078125, 0.234375)] ### 0 is 'lightsalmon', 0.05 is 'crimson'
    elif st == 'TL':
        colors_1 = [(0.390625, 0.58203125, 0.92578125), (0, 0, 0.5)] ### 0 is 'lightskyblue', 0.05 is 'navy'
    else:
        colors_1 = [(0.578125, 0, 0.82421875), (0.86328125, 0.625, 0.86328125)] ### 0 is 'darkviolet', 0.05 is 'plum'
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    if colors_2 == 'None': colors_2 = [(0.75, 0.75, 0.75), (1, 1, 1)]
    # create the second colorscale called "grey_scale"
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N *(1-cutoff)/cutoff)
    #cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N * 9)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)]
    cmap_list += [cmap_2(i) for i in range(cmap_2.N)] 
    cmap_bounds = np.append(np.linspace(0,cutoff,100), np.linspace(cutoff,1,100*int((1-cutoff)/cutoff)))
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, yticklabels = ylabels) 
    # without linewidth=.5, the figure looks much better for section 1
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))


def regular_heatmap(matrix, out_path, out_fname, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0, vmin=0, vmax = 10):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    f, ax = plt.subplots(figsize=figsize)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(0.99609375, 0.625, 0.4765625), (0.859375, 0.078125, 0.234375)] ### 0 is 'lightsalmon', 0.05 is 'crimson'
    elif st == 'TL':
        colors_1 = [(0.390625, 0.58203125, 0.92578125), (0, 0, 0.5)] ### 0 is 'lightskyblue', 0.05 is 'navy'
    else:
        colors_1 = [(0.86328125, 0.625, 0.86328125), (0.578125, 0, 0.82421875)] ### 0 is 'darkviolet', 0.05 is 'plum'
    # create the second colorscale called "grey_scale"
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)] 
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    sns.set(style="white")
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, \
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, yticklabels = ylabels) 
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    plt.close()


def heatmap_hybrid(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0, row_cluster=True, col_cluster=False):
    cluster_mask = matrix.isnull()
    cluster_map = sns.clustermap(matrix.fillna(0), mask = cluster_mask, visible=False, row_cluster=row_cluster, col_cluster=col_cluster)
    row_order = cluster_map.dendrogram_row.reordered_ind
    new_index = [matrix.index[i] for i in row_order]
    matrix_reorder = matrix.reindex(new_index)
    sns.set(font_scale = 2)
    two_scale_heatmap(matrix_reorder, out_path, out_fname, cutoff = cutoff, mask=mask, st = st, \
                      colors_1 = colors_1, colors_2 = colors_2, ylabels = ylabels, \
                      map_square = map_square, figsize = figsize, yaxis_rot = yaxis_rot)  


def metric_col(st, age1, age2):
    if st == 'Prop':
        fc_coln = 'Propensity_%s-%s_FC'%(age1, age2)
        pval_coln = 'Propensity_%s-%s_pval'%(age1, age2)
        st_filter = 'AGG'
    else:
        fc_coln = '%sv%s_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_pval'%(age1[0], age2[0])
        st_filter = st
    return st_filter, fc_coln, pval_coln


def list_to_val(protein_list, df, st, age1, age2, cutoff = 0.05, key = 'N. furzeri Protein Id', out_coln = None, \
    tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis'], sign = 'both'):
    if st == 'Prop':
        if out_coln == None: fc_coln = '%sv%s_prop_logFC'%(age2[0],age1[0])#'Propensity_%s-%s_FC'%(age1, age2)
        pval_coln = '%sv%s_prop_pval'%(age2[0],age1[0])
        df_out = df[df['Sample.Type']=='AGG']
    else:
        if out_coln == None: fc_coln = '%sv%s_logFC'%(age2[0],age1[0])#'%sv%s_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_pval'%(age2[0], age1[0])
        df_out= df[df['Sample.Type']==st]
    if out_coln != None: fc_coln = out_coln
    if sign.lower() == 'pos':
        df_out = df_out[(df_out[pval_coln]<=cutoff) & (df_out[fc_coln]>0)]
    elif sign.lower() == 'neg':   
        df_out = df_out[(df_out[pval_coln]<=cutoff) & (df_out[fc_coln]<0)]
    else:
        df_out = df_out[df_out[pval_coln]<=cutoff]
    df_out = df_out.pivot(index = key, columns = 'Tissue', values = [fc_coln])
    df_out.columns = df_out.columns.get_level_values(1)
    df_out = df_out.reindex(protein_list)
    return df_out


######## Section 3: Given the protein list in order, obtain the other parameters and plot two_scale_heatmap #######
### based on the sorted significant hits, overlay their logFC to generate heatmap 
z_cutoff = None
df_all = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
fpath = 'Overlap/Heatmap'
if not os.path.exists(fpath): os.makedirs(fpath)
sign = 'Pos' # if 'UP' use vmin=0 and vmax=2,  if 'DOWN', use vmin=-2 and vmax = 0
for age_pair in [('Young','Old'), ('Old','TERT'), ('Young','TERT')]:
    age1, age2 = age_pair
    age_comp = '%sv%s'%(age2[0], age1[0])
    # read the organized heatmap order
    for st_query in ['AGG','TL','Prop']:
        query_fname = '%sSig%sZ%s_%s_pval_sorted.csv'%(st_query, sign, z_cutoff, age_comp) if z_cutoff !=None else '%sSig%s_%s_pval_sorted.csv'%(st_query,sign, age_comp)
        df_query = pd.read_csv(os.path.join(fpath, query_fname))
        protein_list = list(df_query['N. furzeri Protein Id'])
        df_query_val = list_to_val(protein_list, df_all, st_query, age1, age2, sign = sign)
        print (df_query_val.min())
        ## Generate a mask for the upper triangle
        print ('%s has a shape of %s x %s'%(st_query, df_query.shape[0], df_query.shape[1]))
        if sign == 'Pos':
            regular_heatmap(df_query_val, fpath, query_fname.replace('.csv', '_heatmap.pdf'), st = st_query, vmin=0, vmax=2)
        elif sign == 'Neg':
            regular_heatmap(df_query_val, fpath, query_fname.replace('.csv', '_heatmap.pdf'), st = st_query, vmin=-4, vmax=0)        
################################################# END OF SECTION 3 ################################################
