#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 12:26:22 2019

@author: yiwenchen
"""

import pandas as pd
import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
from scipy import stats
from sklearn import svm
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
#from sklearn.decomposition import PCA

mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 4)
sns.set_style("white")


def make_kde(*args, **kwargs):    
    sns.kdeplot(*args, cmap=next(make_kde.cmap_cycle), **kwargs)    
    
def spearmanfunc(x, y, **kwargs):
    rho, pval = stats.spearmanr(x, y)
    ax = plt.gca()
    ax.annotate('r = %.2f, p-value = %.2f'%(rho, pval), xy=(.1, .9), xycoords = ax.transAxes)   

# make pair-wise scatter plot 
def grid_correlomics(fig_name, df_data, matrix_coln, info_coln, hue_coln = None, fig_path ='', hue_order = False, colors_list = False, alpha = 1):
    colns = matrix_coln + info_coln 
    df_plot = df_data[colns].dropna()
    if colors_list == False: 
        colors = sns.color_palette('husl', 10)
        colors_list = colors.as_hex()
    sns.set(style = "ticks")
    if hue_coln == None:
        g = sns.PairGrid(df_plot, palette = sns.color_palette(colors_list))
    else:
        if hue_order == False:
            g = sns.PairGrid(df_plot, hue = hue_coln, vars = matrix_coln, palette = sns.color_palette(colors_list), diag_sharey= False)
        else:
            g = sns.PairGrid(df_plot, hue = hue_coln, vars = matrix_coln, hue_order = hue_order, palette = sns.color_palette(colors_list),diag_sharey= False)
    # plot histogram in the diagonal space
    g = g.map_diag(plt.hist, bins=10, histtype = 'step', linewidth = 3)#, density = True) 
#    for i, coln in enumerate(matrix_coln):
#        hist_samebin(coln, df_data, g.axes[i,i], hue_coln, hue_order = hue_order, colors_list = colors_list, normsum = True)
#        g.axes[i,i].set_ylim((0,1))
    # plot scatter plot on the upper corner
    g = g.map_upper(plt.scatter, alpha = alpha)
    make_kde.cmap_cycle = itertools.cycle([sns.light_palette(i, as_cmap=True) for i in colors_list])
    # plot kde in the lower corner (display spearman coefficient if appropriate)
    g = g.map_lower(make_kde)
    if hue_coln == None:
        g.map_lower(spearmanfunc)
    elif df_data[hue_coln].nunique() == 1:
        g.map_lower(spearmanfunc)
    g = g.add_legend() 
    g.savefig(os.path.join(fig_path, fig_name+'.png'))
    plt.close()

# find the seq given the protein id
def search_fasta(pid, seqs, pnames):
    if pid not in pnames:
        return 'no matched protein found'
    else:
        return seqs[pnames.index(pid)]

# ################## Section 1: Initialize the matrix/table if not available ##################
# seeding_dir = 'SysSeeding'
# seeding_fname = 'SysSeedingAll.csv'
# if os.path.exists(os.path.join(seeding_dir, seeding_fname.replace('.csv', '_metrics.csv'))) == False:
#     df_data = pd.read_csv(os.path.join(seeding_dir, seeding_fname))
#     df_MS = pd.read_csv(os.path.join('','kf_combined_prop.csv'))
#     df_data = df_data.merge(df_MS, on = ['N. furzeri Protein Id','Tissue','Sample.Type','N. furzeri Final Symbol','Human'], how = 'left')
#     # merge all the properties as  well as the raw MS data in
#     df_all = pd.read_csv(os.path.join('Properties', 'Killifish_DISOLLRCider.csv'))
#     df_data = df_data.merge(df_all.iloc[:,2:], on=['N. furzeri Protein Id'], how = 'left')
#     df_data.to_csv(os.path.join(seeding_dir, seeding_fname.replace('.csv', '_metrics.csv')),index=False)
# else:
#     df_data = pd.read_csv(os.path.join(seeding_dir, seeding_fname.replace('.csv', '_metrics.csv')))
# #################################### End of Section 1 ####################################


################ Section 1.1: A twist of section 1 to generate charge/hydropathy plot on protein of interest (use aggreation status as hue) ###############
# seeding_dir = 'SysSeeding'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# # plot scatter plot of various physical parameters
# # specify aggregation status color
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na',
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# var_x = 'MeanNetCharge'; var_y = 'uversky_hydropathy'#'Delta_Propensity_Old-Young', log2_Y, MeanNetCharge, YvO_logFC
# df_data.dropna(subset=[var_x,var_y],inplace=True)
# df_data = df_data[df_data['Aggregate']!='TBD']
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# # add the line to separate IDP from globular protein
# ax.set_xlim(-0.01,0.2)
# ax.set_ylim(0.28,0.55)
# x_separator = np.arange(ax.get_xlim()[0],ax.get_xlim()[1]+0.01,0.01)
# ax.plot(x_separator,(x_separator+1.151)/2.785,'k') # check reference at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373528/
# ax.fill_between(x_separator,(x_separator+1.151)/2.785,ax.get_ylim()[1], facecolor='lightblue',alpha=0.5, label='folded')
# ax.fill_between(x_separator,ax.get_ylim()[0],(x_separator+1.151)/2.785, facecolor='lightyellow',alpha=0.5, label='natively unfolded')
# # add the acutal in vivo aggregation data
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_data, ax = ax, \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# #for line in df_data.index:
# ##    x_offset = 0.00; y_offset =-0.001 # this works for #'kappa' as y
# #    x_offset = 0.0001; y_offset =+0.001 # this works for #'Delta_Propensity_Old-Young'
# #    ax.text(df_data[var_x][line]+x_offset, df_data[var_y][line]+y_offset, df_data['Human'][line], horizontalalignment='left', size='small', color='black', weight='regular')
# ax.legend(markerscale=6, bbox_to_anchor=(1.05, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(seeding_dir, 'SysSeeding_uverskyplot.pdf'))
# ################################################## End of Section 1 ###############################################

# ############### Section 1.2: A twist of section 1 to generate charge/hydropathy plot on protein of interest (use aggreation status as hue) ###############
# seeding_dir = 'SysSeeding'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# # plot scatter plot of various physical parameters
# # specify aggregation status color
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na',
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# var_y = 'MeanNetCharge'; var_x = 'uversky_hydropathy'#'Delta_Propensity_Old-Young', log2_Y, MeanNetCharge, YvO_logFC
# df_data.dropna(subset=[var_x,var_y],inplace=True)
# df_data = df_data[df_data['Aggregate']!='TBD']
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# # add the line to separate IDP from globular protein
# ax.set_ylim(-0.01,0.2)
# ax.set_xlim(0.28,0.55)
# x_separator = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 0.001)
# # <H>_boundary = (<R>+1.151)/2.785, or equivalently <R> = 2.785<H>-1.151
# ax.plot(x_separator,(2.785*x_separator-1.151),'k') # check reference at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373528/
# ax.fill_between(x_separator,(2.785*x_separator-1.151),ax.get_ylim()[1], facecolor='lightblue',alpha=0.5, label='folded')
# ax.fill_between(x_separator,ax.get_ylim()[0],(2.785*x_separator-1.151), facecolor='lightyellow',alpha=0.5, label='natively unfolded')
# # add the acutal in vivo aggregation data
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_data, ax = ax, \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# #for line in df_data.index:
# ##    x_offset = 0.00; y_offset =-0.001 # this works for #'kappa' as y
# #    x_offset = 0.0001; y_offset =+0.001 # this works for #'Delta_Propensity_Old-Young'
# #    ax.text(df_data[var_x][line]+x_offset, df_data[var_y][line]+y_offset, df_data['Human'][line], horizontalalignment='left', size='small', color='black', weight='regular')
# ax.legend(markerscale=6, bbox_to_anchor=(1.05, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(seeding_dir, 'SysSeeding_uverskyplot_flip.pdf'))
# ################################################# End of Section 1 ###############################################


##################### Section 2: Output the fasta sequence on protein of interests ####################
# seeding_dir = 'SysSeeding'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# sys_protein = df_data['N. furzeri Protein Id'].tolist()
# # read in the fasta file
# fasta_dir = 'Properties'
# pfasta_table = 'longest-Protein-seqs_Nfu-ncbi100_andMito.fasta'
# f_out = open(os.path.join(seeding_dir, 'SysSeeding.fasta'),'w')
# f = open(os.path.join(fasta_dir, pfasta_table), 'r')
# # seq files
# seqs = f.read().split('\n>')
# #pnames = [i.split('\n')[0].split('|')[2] for i in seqs if '|' in i else i.split('\n')[0].split(' ')[0]]
# pnames = [i.split('\n')[0].split('|')[2] if '|' in i else i.split('\n')[0].split(' ')[0] for i in seqs]
# for pid in sys_protein:
#     out = search_fasta(pid, seqs, pnames)
#     f_out.write('%s\n'%out)
# f_out.close()
########################################### End of Section 4 #############################################

    
################# Section 3: Plot scatterplot on protein of interest (use aggreation status as hue) ###############
# # plot scatter plot of various physical parameters
# # specify aggregation status color
# seeding_dir = 'SysSeeding'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na',
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# var_x = 'OvY_prop_logFC'; var_y = 'delta'#, log2_Y, MeanNetCharge, YvO_logFC
# df_data = df_data[df_data['Aggregate']!='TBD']
# df_data.dropna(subset=[var_x,var_y],inplace=True)
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_data, ax = ax, \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# for line in range(0,df_data.shape[0]):
#     x_offset = 0.00; y_offset =-0.001 # this works for #'kappa' as y
# #    x_offset = 0.0001; y_offset =+0.02 # this works for #'Delta_Propensity_Old-Young'
# #    ax.text(df_data[var_x][line]+x_offset, df_data[var_y][line]+y_offset, df_data['Human'][line], horizontalalignment='left', size='small', color='black', weight='regular')
# ax.legend(markerscale=6, bbox_to_anchor=(1.4, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(seeding_dir, 'SysSeeding_%sv%s_noannotation.pdf'%(var_y, var_x)))
################################################### End of Section 7 ###############################################

#################### Section 4: Visualize fraction of cells with foci in a ranked dot/bar graphh ##################
# ## read matrix of interest
# seeding_dir = 'SysSeeding'
# seeding_fname = 'SysSeedingAll.csv'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# ##################################### End of Section 1 ####################################
# ## plot scatter plot of various physical parameters
# var_x = 'FracWithFoci'; var_y = 'Human'#'Delta_Propensity_Old-Young', log2_Y, MeanNetCharge, YvO_logFC
# df_data.dropna(subset=[var_x,var_y],inplace=True)
# df_plot = df_data.sort_values(by=[var_x, 'Tissue'], ascending=True)
# ### option 1: make a bar plot
# # Make a fake dataset:
# height = df_plot[var_x]
# bars = df_plot[var_y]
# y_pos = np.arange(len(bars))
# # Create bars
# fig, ax = plt.subplots(1,1, figsize=(20,10))
# ax.bar(y_pos, height)
# ax.set_ylabel('fraction of cells with foci')
# # Create names on the x-axis
# plt.xticks(y_pos, bars, fontsize=20, rotation=90)
# # Show graphic
# fig.savefig(os.path.join(seeding_dir, 'SysSeedingFracWithFoci_bar.pdf'))
# ### option 2: make a dot plot
# # Make a fake dataset:
# height = df_plot[var_x]
# bars = df_plot[var_y]
# y_pos = np.arange(len(bars))
# # Create bars
# fig, ax = plt.subplots(1,1, figsize=(20,10))
# ax.scatter(y_pos, height, s=500)
# ax.set_ylabel('fraction of cells with foci')
# # Create names on the x-axis
# plt.xticks(y_pos, bars, fontsize=20, rotation=90)
# # Show graphic
# fig.savefig(os.path.join(seeding_dir, 'SysSeedingFracWithFoci_dot.pdf'))
# ##### option 3: make a dot plot with dot color indicative of tissue origin
# df_plot = pd.read_csv('SysSeeding/SysSeedingAll_dotplot.csv')
# df_plot.dropna(subset=[var_x,var_y],inplace=True)
# df_plot = df_plot.sort_values(by=[var_x, 'Tissue'], ascending=True)
# height = df_plot[var_x]
# bars = df_plot[var_y]
# c_dict = {'Brain':'#d53e4f', 'Gut':'#fc8d59', 'Heart':'#fee08b',\
#                   'Liver':'#ffffbf', 'Muscle':'#e6f598', 'Skin':'#99d594', 'Testis':'#3288bd'}
# color = [c_dict[i] for i in df_plot['Tissue']]
# y_pos = np.arange(len(bars))
# # Create bars
# fig, ax = plt.subplots(1,1, figsize=(14,4))
# for x, y, c in zip (y_pos, height, color):
#     ax.scatter(x, y, s=200, color =c, edgecolors='k')
# ax.set_ylabel('fraction of cells with foci', fontsize=18)
# plt.yticks(np.arange(0,1.1,0.2), fontsize=18, rotation = 90)
# # Create names on the x-axis
# plt.xticks(y_pos, bars, fontsize=18, rotation=90)
# # Show graphic
# fig.savefig(os.path.join(seeding_dir, 'SysSeedingFracWithFoci_dot_tissue.pdf'))
############################
