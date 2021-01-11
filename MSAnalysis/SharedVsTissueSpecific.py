#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:11:46 2019

this file is used to plot a grid matrix of heatmap where only the lower half is shown

@author: yiwenchen
"""
from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import itertools
import matplotlib as mpl

#import ProteinBarChart_YCCT as kfplt

sns.set(style="white")
matplotlib.rcParams['pdf.fonttype'] = 42

# specify a function that compare pair of list and output the number of overlap in each 
def overlap_matrix(var_list, var_dict):
    matrix = np.zeros(shape = (len(var_list), len(var_list)), dtype=int)
    for x, y in itertools.combinations(var_list, r=2):   
        x_index = var_list.index(x); y_index = var_list.index(y)
        matrix[x_index, y_index] = len(set(var_dict[x]) & set(var_dict[y]))
        matrix[y_index, x_index] = len(set(var_dict[x]) & set(var_dict[y]))
    for i, var in enumerate(var_list):
        matrix[i,i] = len(set(var_dict[var]))
    df_matrix = pd.DataFrame(data = matrix, index = var_list, columns = var_list)
    return matrix, df_matrix


def half_heatmap(matrix, out_path, out_fname, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = True, \
                      figsize = (11,9), yaxis_rot = 0, vmin=0, vmax = 150):
    if mask == 'None':
        mask = np.zeros_like(matrix, dtype=np.bool)
        mask[np.triu_indices_from(mask, k = 1)] = True
    f, ax = plt.subplots(figsize=figsize)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(1, 1, 1), (0.859375, 0.078125, 0.234375)] ### 0 is 'lightsalmon', vmax is 'crimson'
    elif st == 'TL':
        colors_1 = [(1, 1, 1), (0, 0, 0.5)] ### 0 is 'lightskyblue', vmax is 'navy'
    else:
        colors_1 = [(1, 1, 1), (0.578125, 0, 0.82421875)] ### 0 is 'darkviolet', vmax is 'plum'
    # create the second colorscale called "grey_scale"
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)] 
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('count_cmap', cmap_list, cmap_1.N)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, linewidths=.5, \
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, yticklabels = ylabels, annot=True, fmt = 'd') 
    ax.set_facecolor('white')
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    plt.close()


def pval_rank(x, delimiter = ',', tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']):
    if delimiter in x:
        out = 0
        for i, tissue in enumerate(tissues):
                if tissue in x: out += 1*(10**i)
    else:
        out = tissues.index(x)
    return out

# this script will be used to show overlap of proteins among different categories

# ####### SECTION 1: FIND THE OVERLAP IN IDENTIFIED PROTEINS AMONG DIFFERENT TISSUES ########
# start with TL or AGG data
fpath = ""
fpathout = "Overlap"
if not os.path.exists(fpathout): os.makedirs(fpathout)
tissues = ["Brain", "Gut", "Heart", "Liver", "Muscle", "Skin", "Testis"]
tissue_colors = ["#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594", "#3288bd"]
tl_proteins = dict()
agg_proteins = dict()
sigtlpos_proteins = dict()
sigaggpos_proteins = dict()
sigproppos_proteins = dict()
fdata = "kf_combined_prop.csv"
df_table = pd.read_csv(os.path.join(fpath, fdata))
#### part 1: find the protein identified in each tissue and assign them to tissue_proteins
for i, tissue in enumerate(tissues):
  tl_proteins[tissue] = df_table["N. furzeri Protein Id"][(df_table["Tissue"] == tissue) & (df_table["Sample.Type"] == "TL")].unique()
  agg_proteins[tissue] = df_table["N. furzeri Protein Id"][(df_table["Tissue"] == tissue) & (df_table["Sample.Type"] == "AGG")].unique()
  sigtlpos_proteins[tissue] = df_table["N. furzeri Protein Id"][(df_table["Tissue"] == tissue) & (df_table["Sample.Type"] == "TL") & (df_table["OvY_logFC"]>=0) & (df_table["OvY_pval"] <=0.05)]
  sigaggpos_proteins[tissue] = df_table["N. furzeri Protein Id"][(df_table["Tissue"] == tissue) & (df_table["Sample.Type"] == "AGG") & (df_table["OvY_logFC"]>=0) & (df_table["OvY_pval"] <=0.05)]
  sigproppos_proteins[tissue] = df_table["N. furzeri Protein Id"][(df_table["Tissue"] == tissue) & (df_table["Sample.Type"] == "AGG") & (df_table["OvY_prop_logFC"]>= 0) & (df_table["OvY_prop_pval"]<=0.05)]
fout = os.path.join(fpathout, fdata.replace("_Combined_prop.csv", "_venn7.pdf"))
#### part 2: generate the matrix
foutpath = 'Overlap'
tl_overlap, df_tl_overlap = overlap_matrix(tissues, tl_proteins)
agg_overlap, df_agg_overlap = overlap_matrix(tissues, agg_proteins)
for st, matrix in zip(['TL','AGG'],[df_tl_overlap, df_agg_overlap]):
    out_fname = '%s_OvY_Overlap.pdf'%st
    half_heatmap(matrix, foutpath, out_fname, st = st, vmax=3600)
# overlap among significant hits (those where OvY_logFC>0)
sigtlpos_overlap, df_sigtlpos_overlap = overlap_matrix(tissues, sigtlpos_proteins)
sigaggpos_overlap, df_sigaggpos_overlap = overlap_matrix(tissues, sigaggpos_proteins)
sigproppos_overlap, df_sigproppos_overlap = overlap_matrix(tissues, sigproppos_proteins)
for st, matrix in zip(['TL','AGG','Prop'],[df_sigtlpos_overlap, df_sigaggpos_overlap, df_sigproppos_overlap]):
    out_fname = '%s_OvY_sigpos_Overlap.pdf'%st
    half_heatmap(matrix, foutpath, out_fname, st = st)
# #################################### END OF SECTION 1 #####################################


# ########### SECTION 2: CREATE MY OWN SORTED LIST FOR TL/AGG/Prop hits by tissue-specificity (Param's output p-val is wrong)########
fpath = ""
#fpathout = "Overlap"
#if not os.path.exists(foutpath): os.makedirs(foutpath)
fdata = "kf_combined_prop.csv"
z_cutoff = None # None or 0.75
fpathout = 'Overlap/Heatmap'
if not os.path.exists(fpathout): os.makedirs(fpathout)
sign = 'Pos'# 'Pos' or 'Neg'
for age_pair in [('Young','Old'), ('Old','TERT'), ('Young','TERT')]:
    age1, age2 = age_pair
    age_comp = '%sv%s'%(age2[0], age1[0])
    key_list = ['N. furzeri Protein Id','N. furzeri Final Symbol','Tissue','Sample.Type','Human']
    tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    for metric in ['TL','AGG','Prop']: #['AGG','Prop']:
        if metric == 'AGG':
            add_on = ['%s_logFC'%age_comp, '%s_pval'%age_comp, '%s_logFC_zscore'%age_comp]
        elif metric == 'Prop':
            add_on = ['%s_prop_logFC'%age_comp,'%s_prop_pval'%age_comp,'%s_prop_zscore'%age_comp]
        st = 'TL' if metric == 'TL' else 'AGG'
        fname_hit = '%sSigZ%s_%s.csv'%(metric, z_cutoff, age_comp) if z_cutoff !=None else '%sSig_%s.csv'%(metric, age_comp)
        type_hit = '%sSig%sZ%s_%s'%(metric, sign, z_cutoff, age_comp) if z_cutoff !=None else '%sSig%s_%s'%(metric, sign, age_comp)
        df_hit = pd.read_csv(os.path.join('Correlomics/Comp', fname_hit))
        df_hit = df_hit[df_hit['Hit Type']==type_hit]
        df_hit_byprot = df_hit.groupby(['N. furzeri Protein Id','N. furzeri Final Symbol'])['Tissue'].apply(','.join).reset_index()
        df_hit_byprot['pval_rank'] = df_hit_byprot['Tissue'].apply(pval_rank)
        df_hit_byprot = df_hit_byprot.sort_values(by='pval_rank', ascending=True)
        df_hit_byprot.to_csv(os.path.join(fpathout, '%s_pval_sorted.csv'%type_hit),index=False)


# ########## SECTION 3: FIND THE OVERLAP IN IDENTIFIED PROTEINS AMONG DIFFERENT TISSUES ########
## start with TL or AGG data
z_cutoff = None
age1 = 'Young'; age2 = 'Old'; age_comp = '%sv%s'%(age2[0], age1[0])
fpathout = 'Overlap'
if not os.path.exists(fpathout): os.makedirs(fpathout)
key_list = ['N. furzeri Protein Id','N. furzeri Final Symbol','Tissue','Sample.Type','Human']
key_list += ['log2_%s_TL'%age for age in [age1, age2]]
key_list += ['log2_%s_zscore_TL'%age for age in [age1, age2]]
tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
tissue_colors = {'Brain':'#d53e4f', 'Gut':'#fc8d59', 'Heart':'#fee08b', 'Liver':'#ffffbf', 'Muscle':'#e6f598', 'Skin':'#99d594', 'Testis':'#3288bd'}
df_raw = pd.read_csv(os.path.join('','kf_combined_prop.csv'))
df_tl = df_raw[df_raw['Sample.Type']=='TL']
df_agg = df_raw[df_raw['Sample.Type']=='AGG']
df_tl_multiple = df_tl['N. furzeri Protein Id'].value_counts()[df_tl['N. furzeri Protein Id'].value_counts()>1] # this list all proteins that have been detected in more than 1 tissues
df_agg_multiple = df_agg['N. furzeri Protein Id'].value_counts()[df_agg['N. furzeri Protein Id'].value_counts()>1] # this list all proteins that have been detected in more than 1 tissues
for metric in ['AGG','Prop']:
    if metric == 'AGG':
        add_on = ['%sv%s_logFC'%(age2[0], age1[0]), '%sv%s_pval'%(age2[0], age1[0]), '%sv%s_logFC_zscore'%(age2[0], age1[0])]
    elif metric == 'Prop':
        add_on = ['%sv%s_prop_logFC'%(age2[0], age1[0]),'%sv%s_prop_pval'%(age2[0], age1[0]),'%sv%s_prop_logFC_zscore'%(age2[0], age1[0])]
#        add_on = ['Delta_Propensity_%s-%s'%(age2, age1),'Propensity_%s-%s_pval'%(age2, age1),'Delta_Propensity_%s-%s_zscore'%(age2, age1)]
    coln_list = key_list + add_on
    fname_hit = '%sSigZ%s_%s.csv'%(metric, z_cutoff, age_comp) if z_cutoff !=None else '%sSig_%s.csv'%(metric, age_comp)
    type_hit = '%sSigPosZ%s_%s'%(metric, z_cutoff, age_comp) if z_cutoff !=None else '%sSigPos_%s'%(metric, age_comp)
    df_hit = pd.read_csv(os.path.join('Correlomics/Comp', fname_hit))
    df_hit = df_hit[df_hit['Hit Type']==type_hit]
    df_hit['N. furzeri Protein Id'].value_counts().to_csv(os.path.join(fpathout,'%s_InMultipleByProtein.csv'%type_hit))
    print (pd.pivot_table(df_hit[['N. furzeri Protein Id','Tissue']],index=['Tissue'],aggfunc='count'))
    # list of hits that is also expressed in mutliple tissues (tl and agg value available)
    hit_multiple = set(df_hit['N. furzeri Protein Id'].values) & set(df_tl_multiple.index) & set(df_agg_multiple.index)
    print (pd.pivot_table(df_hit[['N. furzeri Protein Id','Tissue']][~df_hit['N. furzeri Protein Id'].isin(hit_multiple)],index=['Tissue'],aggfunc='count'))
    df_hit_multiple = df_agg[coln_list][df_agg['N. furzeri Protein Id'].isin(hit_multiple)]
    # this dataframe is used to track expression differences across tissues
    df_hit_multiple_expY = pd.pivot_table(df_hit_multiple[coln_list], index = ['N. furzeri Protein Id','N. furzeri Final Symbol'], values = ['log2_%s_TL'%age1, add_on[0]], aggfunc=[np.std,np.mean,'count'])
    df_stats = df_hit_multiple_expY['std']/df_hit_multiple_expY['mean']
    fig_tl, ax_tl = plt.subplots(1,3, figsize=(12,4))
    ax_tl[0].hist(df_stats['log2_%s_TL'%age1].dropna().values, bins=50,range=(0,1))
    ax_tl[0].set_title('log2_%s_TL: std/mean'%age1)
    ax_tl[1].hist(df_stats['log2_%s_TL'%age1].dropna().values, bins=50,range=(0,1),cumulative=True, density=True, histtype='step')
    ax_tl[1].set_yticks(np.arange(0,1.1,0.1))
    ax_tl[1].set_title('log2_%s_TL: std/mean'%age1)
    ax_tl[2].hist(df_hit_multiple_expY['count']['log2_%s_TL'%age1].dropna().values, bins = range(1,8))
    ax_tl[2].set_title('# of tissues expressed')
    fig_tl.suptitle('%s_InMultipleTissues'%type_hit)
    fig_tl.savefig(os.path.join(fpathout, '%s_InMultipleTissues.pdf'%type_hit))
    df_hit_mulitple_pivot = pd.pivot_table(df_hit_multiple[coln_list], index = ['N. furzeri Protein Id','N. furzeri Final Symbol'], columns = ['Tissue'])
    # note that the log2_Young_TL is filtered and will only have notna value if the protein is also detected in AGG (I used df_agg as filter)
    df_hit_mulitple_pivot.to_csv(os.path.join(fpathout, '%s_InMultipleTissues.csv'%type_hit))
# ##################################### END OF SECTION 3 ####################################
