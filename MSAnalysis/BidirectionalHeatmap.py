#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 21:55:59 2019

this is used to generate heatmap with bi-rectional color scale (usually either just log2FC or a composite metric of log2FC*-logPval)

@author: yiwenchen
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import numpy as np


mpl.rcParams['pdf.fonttype'] = 42

sns.set(style="white")


## note this pval_filter will only keep term where the tissue has a significant enrichment
def pval_filter(df, tissues, pvalcutoff = 0.05, key = 'Description', direction = 'both'):
    df_new = df.copy(deep=True)
    df_new['sig_tissue'] = ''
    for i, tissue in enumerate(tissues):
        row_reset_i = df[df['pvalue_%s'%tissue]>pvalcutoff].index.tolist() 
        df_new.loc[row_reset_i, [j for j in df_new.columns if tissue in j]] = np.nan
        if direction == 'both':
            df_new.loc[df[df['pvalue_%s'%tissue]<=pvalcutoff].index.tolist(), 'sig_tissue'] += str(i+1)
        elif direction == 'up':
            df_new.loc[df[(df['pvalue_%s'%tissue]<=pvalcutoff) & (df['NES_%s'%tissue]>0)].index.tolist(), 'sig_tissue'] += str(i+1)
        elif direction == 'down':
            df_new.loc[df[(df['pvalue_%s'%tissue]<=pvalcutoff) & (df['NES_%s'%tissue]<0)].index.tolist(), 'sig_tissue'] += str(i+1)
    df_new['sig_tissue'] = pd.to_numeric(df_new['sig_tissue'],errors='coerce')
    df_new.sort_values(by=['sig_tissue'], inplace = True)
    df_new = df_new.drop_duplicates(subset=key)
    return df_new


## note this FC_filter will only keep term where the enrichment is either positive or negative
def FC_filter(df, tissues, fcdir, key = 'Description'):
    df_new = df.copy(deep=True)
    for i, tissue in enumerate(tissues):
        if fcdir == 'UP':
            row_reset_i = df[df['NES_%s'%tissue]<0].index.tolist()
        elif fcdir == 'DOWN':
            row_reset_i = df[df['NES_%s'%tissue]>0].index.tolist()
        df_new.loc[row_reset_i, [i for i in df_new.columns if tissue in i]] = np.nan
    df_new = df_new.drop_duplicates(subset=key)
    return df_new


def custom_heatmap(matrix, out_path, out_fname, mask_color='Gray', vmin = -5, vmax = 5, figsize=(20,10), yrot=0):
    mask = matrix.isnull() # this is mask nan in the matrix
    colormap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])   
    f, ax = plt.subplots(figsize=figsize)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    matrix_map = sns.heatmap(matrix, cmap = colormap, mask=mask, vmin = vmin, vmax = vmax, \
                             square=True, linewidths=.5, ax = ax, cbar_kws={"shrink": .5}, yticklabels = True)
    ax.set_title(out_fname.replace('.pdf',''))
    matrix_map.set_facecolor(mask_color)
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yrot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    plt.close()


def custom_clustermap(matrix, out_path, out_fname, cutoff = 0.05, mask_color='Gray', vmin = -5, vmax = 5, figsize=(4,10)):
    mask = matrix.isnull()
    colormap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])   
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 2)
    cluster_map = sns.clustermap(matrix.fillna(0), cmap = colormap, mask = mask, \
                                col_cluster = False, square = True,\
                                yticklabels = True, figsize = figsize)
    cluster_map.ax_heatmap.set_facecolor(mask_color)
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
#    matrix_map.tick_params(labelsize=12)
#    matrix_map.tick_params(axis='y', labelrotation=0)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    cluster_map.savefig(os.path.join(out_path, out_fname))
    plt.close()


def heatmap_cluster_hybrid(matrix, out_path, out_fname, mask_color='Gray', row_cluster=False, col_cluster=False, vmin = -5, vmax = 5,figsize=(20,10),yrot=0):
    mask = matrix.isnull()
    cluster_map = sns.clustermap(matrix.fillna(0), mask = mask, row_cluster=row_cluster, col_cluster=col_cluster, visible=False)
    if row_cluster == True:
        row_order = cluster_map.dendrogram_row.reordered_ind
        new_index = [matrix.index[i] for i in row_order]
        matrix_reorder = matrix.reindex(new_index)
    if col_cluster == True:
        col_order = cluster_map.dendrogram_col.reordered_ind
        new_col = [matrix.columns[i] for i in col_order]
        if row_cluster != True:
            matrix_reorder = matrix[new_col]
        elif row_cluster == True:
            matrix_reorder = matrix_reorder[new_col]
    if (row_cluster != True) and (col_cluster != True):
        matrix_reorder = matrix.copy(deep=True)
    sns.set(font_scale = 3)
    custom_heatmap(matrix_reorder, out_path, out_fname, mask_color = mask_color, vmin = vmin, vmax = vmax, figsize=figsize, yrot=yrot)
    plt.close()     


################################################### SECTION 7 ################################################
# Simple script where you generate bi-directional heatmap on specific table of interest
# this is applying the p-value cutoff filter to mask out insignificant terms    
age_comp = 'TvO'# can either be 'OvY', 'TvO', 'TvY'
fpath = "GSEA/Cleanedup/%s"%age_comp
foutpath = "GSEA/Cleanedup/Fig2"
if not os.path.exists(foutpath): os.makedirs(foutpath)

metric = 'NES'
st = 'TL'# TL, AGG, Prop
pval_cutoff = 0.05
nbest = 3 # number of term that are either significantly enrichment in top (increase with age) or bottom (decrease with age)
nsig = 0 # minimum number of term that shared significant enrichment to be visualized
#direction = 'both'
directions = ['down', 'up']#,'down', 'both']
for direction in directions:
#    for dtype in ['GOGSEA', 'KEGG', 'KEGG-Modules','MSigDb-GO-ALL','MSigDb-GOCC', 'MSigDb-HALLMARKS', 'DO']:
    for dtype in ['KEGG', 'KEGG-Modules', 'MSigDb-HALLMARKS']:
        fname = 'All_%s_%s_%s_nbest%s_nsig%s_%s_fromraw.csv'%(st, dtype, age_comp, nbest, nsig, direction.lower())
        tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
        if 'T' in age_comp: tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
        df_i = pval_filter(pd.read_csv(os.path.join(fpath, fname)), tissues, pvalcutoff = pval_cutoff, direction = direction)
        #extract the relevant columns
        matrix_coln = ['Description'] + ['%s_%s'%(metric, i) for i in tissues]
        df_map_i = df_i[matrix_coln]
        df_map_i.set_index('Description', inplace = True)
        df_map_i = df_map_i.dropna(how='all')
        if df_map_i.empty == False:
            custom_heatmap(df_map_i, foutpath, '%s_heatmap.pdf'%(fname.replace('.csv', '')))
    #        if df_map_i.shape[0]>1:
    #            heatmap_cluster_hybrid(df_map_i, foutpath,'%s_clustermap.pdf'%(fname.replace('.csv','')),figsize=(60,20), yrot=0,row_cluster=True,col_cluster=False)
################################################### End of SECTION 7 ################################################