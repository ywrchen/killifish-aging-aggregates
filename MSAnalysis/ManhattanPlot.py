#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:45:11 2019

# create an Manhattan style plot to show FC and p-value, vary the y to try out different layout
and consider only plahying with point size etc to illustrate other points
# the essential here is to create scatter plot 

@author: yiwenchen
"""


import pandas as pd
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
  

mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 3)
sns.axes_style("white")

# create an manhattan plot on variable of interest
# x is the variable to be used to map the category, y is the variable on the y-axis, typically -log10pval, 
# s is the variable where the point size is determined from
def manhattan(df, x, y, title, s, ax = '', y_colormap = '', s_colormap = 'default', fpath = '', sizes='', \
              cutoff='None', cutoff_c = '#D4D4D4', legend = 'brief', s_max = 'None', tissues = '', annotate = True):
    sns.set_style("white", {'axes.linewidth':0.5})
    if ax == '':
        fig, ax_m = plt.subplots(1,1, figsize=(20, 10))
    else:
        ax_m = ax
    df= df.sort_values(x)
    if tissues !='': df= df[df['Tissue'].isin(tissues)]
    df['x-axis'] = range(len(df))
    if 'zscore' not in y: ax_m.set_yscale('log', basey=2)
    if cutoff == 'None':
        n_tissue = len(df['Tissue'].unique())
        if sizes =='':
            sns.scatterplot(x='x-axis', y=y, data = df, hue = 'Tissue', size=s, \
                    palette = sns.color_palette(y_colormap[:n_tissue]), ax=ax_m, edgecolor='k') 
        else:
            sns.scatterplot(x='x-axis', y=y, data = df, hue = 'Tissue', size=s, sizes = sizes,\
                    palette = sns.color_palette(y_colormap[:n_tissue]), ax=ax_m, edgecolor='k')
    else:
        cutoff_coln = s
        cutoff_val = cutoff
        df_grey = df[df[cutoff_coln]<cutoff_val]
        df_color = df[df[cutoff_coln]>=cutoff_val]
        n_grey = len(df_grey['Tissue'].unique())
        n_color = len(df_color['Tissue'].unique())
        if s_max == 'None': s_max = df[s].max()
        sns.scatterplot(x='x-axis', y=y, data = df_grey, size=s, hue = 'Tissue', sizes = sizes,\
                    palette = sns.color_palette([cutoff_c]*n_grey), ax=ax_m, edgecolor='lightgray', size_norm = (1,s_max), legend = legend)
        sns.scatterplot(x='x-axis', y=y, data = df_color, size=s, hue = 'Tissue', sizes = sizes,\
                    palette = sns.color_palette(y_colormap[:n_color]), ax=ax_m, edgecolor='k', size_norm = (1,s_max), legend = False)
    if annotate == True:
        for index_r, row in df_color.iterrows():
            ax_m.text(row['x-axis']+0.2, row[y], row['Human'], horizontalalignment='left', size=8, color='black', weight='light')
#        for line in range(0,df.shape[0]):
#            ax.text(df['x-axis'].iloc[line]+0.2, df[y].iloc[line], df['Human'].iloc[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
    ax_m.set_title(title.replace('.pdf', ''))
    ax_m.legend(markerscale=1, bbox_to_anchor=(1,1))
    ax_m.get_xaxis().set_visible(False)
    if 'zscore' not in y:
        ax_m.yaxis.set_ticks([2**-4, 2**-2, 2**0, 2**2, 2**4])
        ax_m.set_ylim([2**-4,2**4])
    else:
        ax_m.set_ylim([-5,5])
    plt.rcParams['ytick.left'] = True
    if ax == '':
        fig.savefig(os.path.join(fpath, title), bbox_inches = 'tight')
        return None
    else:
        return ax 
    

# create an manhattan plot on variable of interest
# x is the variable to be used to map the category, y is the variable on the y-axis, typically -log10pval, 
# s is the variable where the point size is determined from

def manhattan_overlay(df, x, y, title, s, ax = '', y_colormap = '', s_colormap = 'default', fpath = '', sizes='', \
              cutoff='None', cutoff_c = '#D4D4D4', legend = 'brief', s_max = 'None', tissues = '', annotate = True):
    sns.set_style("white", {'axes.linewidth':0.5})
    if ax == '':
        fig, ax_m = plt.subplots(1,1, figsize=(20, 10))
    else:
        ax_m = ax
    if tissues !='': df= df[df['Tissue'].isin(tissues)]
    df['x-axis'] = df.groupby([x]).ngroup()
    ax_m.set_yscale('log', basey=2)
    if cutoff == 'None':
        n_tissue = len(df['Tissue'].unique())
        if sizes =='':
            sns.scatterplot(x='x-axis', y=y, data = df, hue = 'Tissue', \
                    palette = sns.color_palette(y_colormap[:n_tissue]), ax=ax_m, edgecolor='k', alpha =0.5) 
        else:
            sns.scatterplot(x='x-axis', y=y, data = df, hue = 'Tissue', size=s, sizes = sizes,\
                    palette = sns.color_palette(y_colormap[:n_tissue]), ax=ax_m, edgecolor='k', alpha =0.5)
    else:
        cutoff_coln = s
        cutoff_val = cutoff
        df_grey = df[df[cutoff_coln]<cutoff_val]
        df_color = df[df[cutoff_coln]>=cutoff_val]
        n_grey = len(df_grey['Tissue'].unique())
        n_color = len(df_color['Tissue'].unique())
        if s_max == 'None': s_max = df[s].max()
        sns.scatterplot(x='x-axis', y=y, data = df_grey, size=s, hue = 'Tissue', sizes = sizes,\
                    palette = sns.color_palette([cutoff_c]*n_grey), ax=ax_m, edgecolor='lightgray', size_norm = (1,s_max), legend = legend)
        sns.scatterplot(x='x-axis', y=y, data = df_color, size=s, hue = 'Tissue', sizes = sizes,\
                    palette = sns.color_palette(y_colormap[:n_color]), ax=ax_m, edgecolor='k', size_norm = (1,s_max), legend = False)
    if annotate == True:
        for index_r, row in df_color.iterrows():
            ax_m.text(row['x-axis']+0.2, row[y], row['Human'], horizontalalignment='left', size=8, color='black', weight='light')
#        for line in range(0,df.shape[0]):
#            ax.text(df['x-axis'].iloc[line]+0.2, df[y].iloc[line], df['Human'].iloc[line], horizontalalignment='left', size='medium', color='black', weight='semibold')
    ax_m.set_title(title.replace('.pdf', ''))
    ax_m.legend(markerscale=1, bbox_to_anchor=(1,1))
    ax_m.get_xaxis().set_visible(False)
#    ax_m.yaxis.set_ticks([2**-4, 2**-2, 2**0, 2**2, 2**4])
    ax_m.set_ylim([2**-4,2**4])
    plt.rcParams['ytick.left'] = True
    if ax == '':
        fig.savefig(os.path.join(fpath, title), bbox_inches = 'tight')
        return None
    else:
        return ax


########### Section 4.1: Do AGG, TL, and Prop on a Specific Age-comparision #########
# this is vertical layout (potrait version)         
fpath = ''
fname = 'kf_combined_prop.csv'
foutpath = 'Tert'
if not os.path.exists(foutpath): os.makedirs(foutpath)
df_all = pd.read_csv(os.path.join(fpath, fname))
# sort by tissue and only keep the aggregate
sts = ['TL', 'AGG', 'Prop']
age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
# define a set of colors for the tissues
tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd'] 
tissues_order = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis'] 
sns.set_style("white", {'axes.linewidth':0.5})
#### part 1: use fc_z score as the y axix ####
for age_pair in age_pairs:
    age1, age2 = age_pair
    fig, ax = plt.subplots(len(sts),1, figsize=(30, 10*len(sts)), sharex = True, sharey=True)
    for st, ax_i, i in zip(sts, ax.flat, range(len(sts))):
        if st == 'Prop':
            df_data = df_all[df_all['Sample.Type']=='AGG'].sort_values('Tissue')
#            fc_coln = 'Delta_Propensity_%s-%s_zscore'%(age1, age2)
#            pval_coln = 'Propensity_%s-%s_pval'%(age1, age2)
            fc_coln = '%sv%s_prop_logFC_zscore'%(age1[0], age2[0])
            pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
            logpval_coln = '-log10pval_%sv%sProp'%(age1[0], age2[0])
        else:
            df_data = df_all[df_all['Sample.Type']==st].sort_values('Tissue')
            fc_coln = '%sv%s_logFC_zscore'%(age1[0], age2[0])
            pval_coln = '%sv%s_pval'%(age1[0], age2[0])
            logpval_coln = '-log10pval_%sv%s'%(age1[0], age2[0])
        df_data[logpval_coln] = -np.log10(df_data[pval_coln])     
        if i == 0:
            ax_i = manhattan(df_data, 'Tissue', fc_coln, '%s%sv%s_Manhattanv2.pdf'%(st, age1[0], age2[0]), logpval_coln, \
              ax = ax_i, y_colormap = tissues_colors, sizes = (10,200), cutoff = -np.log10(0.05), s_max = 6, annotate=False)
        else:
            ax_i = manhattan(df_data, 'Tissue', fc_coln, '%s%sv%s_Manhattanv2.pdf'%(st, age1[0], age2[0]), logpval_coln, \
              ax = ax_i, y_colormap = tissues_colors, sizes = (10,200), cutoff = -np.log10(0.05), legend = False, s_max = 6, annotate=False)
    fig.savefig(os.path.join(foutpath, '%sv%s_zscore_Manhattanv2.pdf'%(age1[0], age2[0])), bbox_inches = 'tight')     
################# End of Section 4.1 ###############


######### Section 7: Get frac of protein that's significantly changing given certain pval and fold change cuttoff #########
# update on 2020.03.30
fpath = ''
fname = 'kf_combined_prop.csv'
foutpath = 'MSPlots'
if not os.path.exists(foutpath): os.makedirs(foutpath)
df_all = pd.read_csv(os.path.join(fpath, fname))
# sort by tissue and only keep the aggregate
for st in ['TL', 'AGG', 'Prop']:
    if st == 'TL':
        df_agg = df_all[df_all['Sample.Type']==st].sort_values('Tissue')
    else:
        df_agg = df_all[df_all['Sample.Type']=='AGG'].sort_values('Tissue')
    age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
    # define a set of colors for the tissues
    tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
    tissues_order = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    sns.set_style("white", {'axes.linewidth':0.5})
    fig, ax = plt.subplots(len(age_pairs),1, figsize=(30, 10*len(age_pairs)), sharex = True, sharey=True)
    for age_pair, ax_i, i in zip(age_pairs, ax.flat, range(len(age_pairs))):
        age1, age2 = age_pair
        fc_coln = '%sv%s_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_pval'%(age1[0], age2[0])
        df_agg['-log10pval_%sv%s'%(age1[0], age2[0])] = -np.log10(df_agg[pval_coln])
        if i == 0:
            ax_i = manhattan(df_agg, 'Tissue', fc_coln, '%s%sv%s_Manhattanv2.pdf'%(st, age1[0], age2[0]), '-log10pval_%sv%s'%(age1[0], age2[0]), \
              ax = ax_i, y_colormap = tissues_colors, fpath = 'PosPropFC', sizes = (10,200), cutoff = -np.log10(0.05), s_max = 6, tissues = tissues_order, annotate = False)
        else:
            ax_i = manhattan(df_agg, 'Tissue', fc_coln, '%s%sv%s_Manhattanv2.pdf'%(st, age1[0], age2[0]), '-log10pval_%sv%s'%(age1[0], age2[0]), \
              ax = ax_i, y_colormap = tissues_colors, fpath = 'PosPropFC', sizes = (10,200), cutoff = -np.log10(0.05), legend = False, s_max = 6,  tissues = tissues_order, annotate = False)
        ax_i.hlines(y=1.5, xmin = 0, xmax = df_agg.shape[0], linewidth=2, color='k', linestyles='dashed')
        ax_i.hlines(y=0.5, xmin = 0, xmax = df_agg.shape[0], linewidth=2, color='k', linestyles='dashed')
    fig.savefig(os.path.join(foutpath, '%s_Manhattanv2_hline.pdf'%st), bbox_inches = 'tight')
############### End of Section 7 ###############
