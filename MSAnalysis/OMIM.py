#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 12:19:12 2019

this is file used to filter map out OMIM entry for proteins with age-assoicated changes
in aggregation propensity

@author: yiwenchen
"""

#import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import math
import scipy.stats
import matplotlib
import seaborn as sns
#import itertools
import matplotlib.backends.backend_pdf
from PyPDF2 import PdfFileWriter, PdfFileReader
import scipy.stats as stats

# import custom functions
#from ProteinBarChart_YCCT import bar_jitter, bar_jitter_all, bar_jitter_panel, bar_ageprop
from ProteinBarChart import bar_sample, bar_MS, bar_MS_separate

sns.set(font_scale = 2)
sns.set_style("white")

# read a file and save the list of protein
def read_list(fname, path, identifier):
    df_protein_list = pd.read_csv(os.path.join(path, fname))
    df_clean = df_protein_list.drop_duplicates()
    p_list = df_clean[identifier].dropna().tolist()
    return p_list


# read a txt file and make each unique gene symbols 
def read_split_csv(fname, split_coln, map_colns, separator = ','):
    df = pd.read_csv(fname, sep='\t', skiprows=3)
    df.dropna(axis=0, subset=[split_coln], inplace=True)
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
    The values in the other columns are duplicated across the newly divided rows.
    '''
    row_accumulator = []
    def splitListToRows(row, separator):
        if separator in row[split_coln]:
            split_row = row[split_coln].split(separator)
        else:
            split_row = [row[split_coln]]
        for s in split_row:
            new_row = row.to_dict()
            new_row[split_coln] = s.strip()
            row_accumulator.append(new_row)
    df.apply(splitListToRows, axis=1, args = (separator, ))
    df_out = pd.DataFrame(row_accumulator)   
    return df_out

# make histogram from dataframe with hue
def hist_hue(df_data, hist_colns, hist_fname, hue_coln, hue_order = None, colors_list = None, hit_count=None, alpha = 1, 
             sharey = False, histtype = 'step', bins= 10, lw = 3, normsum = False, ks_test = False, directory = ''):
    if hue_order == None:
        hue_order = sorted(set(df_data[hue_coln].values))
        if 'null' in hue_order: hue_order.append(hue_order.pop(hue_order.index('null')))
        if 'detected' in hue_order: hue_order.append(hue_order.pop(hue_order.index('detected')))
    if (hue_order != None) & (colors_list != None):
        if len(hue_order) != len(colors_list):
            print ('Default colors used instead because the specified color do not match with hue_coln')
            colors_list = sns.color_palette('husl', len(hue_order))
    if colors_list == None:
        colors_list = sns.color_palette('husl', len(hue_order))
    hist_coln_num = 4;     
    hist_row_num = int(math.ceil(len(hist_colns)/1./hist_coln_num))
    fig_hist, ax_hist = plt.subplots(nrows = hist_row_num, ncols = hist_coln_num, figsize = (16,4*hist_row_num), sharey = sharey)
    ks_test_out = [[],[], []]
    for ax, metric in zip(ax_hist.flat, hist_colns):
        df_plot = df_data[[metric,hue_coln]].dropna(subset=[metric]).apply(pd.to_numeric,errors='coerce')
        metric_min = df_plot[metric].min()
        metric_max = df_plot[metric].max()
        for hue_i, color in zip(hue_order, colors_list):
            hist_x = df_plot[metric][df_data[hue_coln]==hue_i]
            if hit_count == hue_i: n_hit = len(hist_x)
            if normsum == False:
                fig = sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins = bins,\
                         hist_kws = {'histtype':histtype,'alpha':alpha,'lw':lw,'color':color,\
                                     'range':(metric_min, metric_max),'label':hue_i})#,'density': True})
            # normsum is to normalize the sum of the bin counts to total of 1 (i.e. histogram on fraction of counts)
            elif normsum == True:
                weights = np.ones_like(hist_x)*1./len(hist_x)
                fig = sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins=bins,\
                             hist_kws = {'histtype':histtype,'alpha':alpha,'lw':lw,'color':color,\
                                         'label':hue_i, 'weights':weights, 'range':(metric_min, metric_max)})
        if ks_test == True:
            if df_data[hue_coln].nunique() == 2:
                hue_col_vs = df_data[hue_coln].unique().tolist()
                ks_stat, ks_pvalue = stats.ks_2samp(df_data[metric][df_data[hue_coln]==hue_col_vs[0]],df_data[metric][df_data[hue_coln]==hue_col_vs[1]])
#                fig.text(0.5, 0.9, 'ks: %.2f, p-value: %.2f'%(ks_stat, ks_pvalue), horizontalalignment='center', size='medium', color='black', transform=ax.transAxes)
                if hit_count!=None:
                    ax.set_title('%s\n ks = %.2f, p-value = %.2f\nn_POI = %s'%(metric, ks_stat, ks_pvalue, n_hit))
                else:
                    ax.set_title('%s\n ks = %.2f, p-value = %.2f'%(metric, ks_stat, ks_pvalue))
                ks_test_out[0].append(metric)
                ks_test_out[1].append(ks_stat)
                ks_test_out[2].append(ks_pvalue)
        else:
            ax.set_title(metric)
    handles,labels = ax.get_legend_handles_labels()
    fig_hist.legend(handles, labels, loc = 'center right')
    plt.tight_layout()
    fig_hist.savefig(os.path.join(directory,'%s.pdf'%hist_fname))
    # retrieve the ks test results
    return ks_test_out


###### this generates a regular heatmap that go from 0 to 1 in one directional color scale ######
def regular_heatmap(matrix, out_path, out_fname, mask='None', st = 'None', annot=False,\
                      colors_1 = 'None', colors_2 = 'None', xlabels = False, ylabels = False, map_square = True, \
                      figsize = (10,10), yaxis_rot = 0, vmin=0, vmax = 1, ax = 'None'):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    if ax == 'None':
        f, ax = plt.subplots(figsize=figsize)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(1, 1, 1), (0.859375, 0.078125, 0.234375)] ### (0.99609375, 0.625, 0.4765625) is 'lightsalmon', (0.859375, 0.078125, 0.234375) is 'crimson'
    elif st == 'WCL':
        colors_1 = [(1, 1, 1), (0, 0, 0.5)] ### (0.390625, 0.58203125, 0.92578125) is 'lightskyblue', (0, 0, 0.5) is 'navy'
    else:
        colors_1 = [(1, 1, 1), (0.578125, 0, 0.82421875)] ### (0.86328125, 0.625, 0.86328125) is 'darkviolet', (0.578125, 0, 0.82421875) is 'plum'
    # create the second colorscale called "grey_scale"
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)] 
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    # Draw the heatmap with the mask and correct aspect ratio
#    sns.set(font_scale = 3)
    if annot == False:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, \
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, xticklabels = xlabels, yticklabels = ylabels) 
    else:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, annot=True, fmt='.2g',\
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, xticklabels = xlabels, yticklabels = ylabels)
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    if ax == 'None': 
        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    else:
        ax.set_title(out_fname.replace('.pdf',''), fontsize=24)
        return ax

def two_scale_heatmap(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', xlabels = False, ylabels = False, \
                      map_square = False, figsize = (10,10), yaxis_rot = 0, ax = None, annot=False):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    if ax == None: f, ax = plt.subplots(figsize=figsize)
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(0.859375, 0.078125, 0.234375), (0.99609375, 0.625, 0.4765625)] ### 0 is 'lightsalmon', 0.05 is 'crimson'
    elif st == 'WCL':
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
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1, annot=annot,\
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, xticklabels = xlabels, yticklabels = ylabels) 
    # without linewidth=.5, the figure looks much better for section 1
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    if ax == None:
        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    else:
        return ax


################## Section 1: Read phenotype map and melt it ##################
# phenotype_list = 'Download/genemap2.txt'
# input_path = 'OMIM'
# protein_map = read_split_csv(os.path.join(input_path, phenotype_list), 'Gene Symbols', ['Mim Number', 'Phenotypes'])
# protein_map.to_csv(os.path.join(input_path, phenotype_list.replace('.txt', '_melt.csv')), index = False)
############################## End of Section 1 ##############################


#matplotlib.rcParams.update({'font.size': 6})
#age_tag_dict = {'Young': ['Young_X128_C', 'Young_X128_N', 'Young_X129_N'],\
#                'Old': ['Old_X126', 'Old_X127_C', 'Old_X127_N'], \
#                'Tert': ['TERT_X129_C', 'TERT_X130_C', 'TERT_X130_N']}
#young_hex = '#91d6e4'; old_hex = '#488cca'; tert_hex = '#7f7f7f'
#tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis']


##### Section 2: filter out proteins with significant age-associated changes in aggregation propensity #####
# updated on 2020.05.16
age1 = 'Young'; age2 = 'Old'; metric = 'AGG'; direction='Pos' # or 'AGG' or 'WCL'
p_cutoff = 0.05; z_cutoff = 0;
fpath_in = 'Correlomics/Comp'
fname_in = '%sSig_%sv%s.csv'%(metric, age2[0], age1[0]) if z_cutoff == 0 else '%sSigZ%s_%sv%s.csv'%(metric, z_cutoff, age2[0], age1[0])
fname_out = '%sSig%s_%sv%s'%(metric, direction, age2[0], age1[0]) if z_cutoff == 0 else '%sSig%sZ%s_%sv%s'%(metric, direction, z_cutoff, age2[0], age1[0])
fpath_out = 'OMIM/SigHits'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# import the hit table
df_hit = pd.read_csv(os.path.join(fpath_in, fname_in))
df_sig = df_hit[df_hit['Hit Type'].str.contains(direction)].dropna(axis=0, subset=['Human'])
protein_map = pd.read_csv(os.path.join('OMIM', 'genemap2_melt.csv'))
# merge all the entries associated with the same term from different tissue in one
df_sig_agg = df_sig.groupby(['N. furzeri Final Symbol', 'N. furzeri Protein Id', 'Human']).agg({'Tissue': ','.join}).reset_index()
#df_sigprop_merge = df_sigprop.groupby(['N. furzeri Final Symbol', 'Human', prop_pval, prop_FC]).agg({'Tissue': ','.join}).reset_index()
df_merge = df_sig.merge(protein_map, how = 'left', left_on = 'Human', right_on = 'Gene Symbols')
df_merge_agg = df_sig_agg.merge(protein_map, how = 'left', left_on = 'Human', right_on = 'Gene Symbols')
##### part 1. update the merged hit table and OMIM entries
df_merge.to_csv(os.path.join(fpath_out, '%s_OMIM.csv'%fname_out),index=False)
df_merge_agg.to_csv(os.path.join(fpath_out, '%s_OMIM_jointissue.csv'%fname_out),index=False)
#### part 2. plot the raw data on WCL, AGG, and Prop for hits with OMIM entries and showed up in multiple tissues
df_multi = df_merge_agg.dropna(axis=0, subset=['Phenotypes'])
df_multi = df_multi[df_multi['Tissue'].str.contains(',')]
df_raw = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
out_pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(fpath_out, '%s_OMIM_multissue.pdf'%fname_out))
key = 'N. furzeri Protein Id'
fig_list = []
#all_proteins = protein_map['Gene Symbols'].tolist()
for pid, protein_hm, tissues in zip(df_multi[key].values, df_multi['Human'].values, df_multi['Tissue'].values):
    tissues = tissues.split(',')
    char_break = 120*len(tissues)
    disease = str(df_multi['Phenotypes'][df_multi[key]==pid].values[0])
    disease = ''.join([char if i % char_break !=(char_break-1) else char+'\n' for i, char in enumerate(disease)])
    fig, ax = bar_MS_separate(df_raw, pid, key, stype=['WCL','AGG','Prop'], tissues=tissues,\
                              fname_out=disease, fpath_out=fpath_out, ind = False, highlight = tissues, sharey=False)
    fig_list.append(fig)
    out_pdf.savefig(fig)
out_pdf.close()
############################## End of Section 2 ##############################
