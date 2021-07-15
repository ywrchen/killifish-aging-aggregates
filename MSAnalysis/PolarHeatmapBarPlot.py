#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 11:43:07 2020

@author: yiwenchen
"""

# make a custom heatmap style bar plot in polar coordinates
# point is better visaluze tissue-level differences

import pandas as pd
import numpy as np
import os
import scipy.stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
 

mpl.rcParams['pdf.fonttype'] = 42
sns.set(style="white")


############ Section 1: Calculate the extent of changes at increasing z-score cutoff ############
# look through matrix then loop through different age comparisons
# either with or without z_cutoff #
# fpath_out = 'Tert/ProliferationIndex'
# pcutoff = 0.05
# age_pairs = [('Young','Old'), ('Old','Tert')]
# z_cutoff = None#None # 0.75
# df_table = pd.read_csv(os.path.join('', 'kf_Combined_prop.csv'))
# tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# sign = 'both'#'both'
# #### part 1: age-associated changes comparison between two age pairs, updated on 2020.06.01
# for metric in ['TL', 'AGG', 'Prop']:
#     st = metric
#     if metric == 'Prop':
#         pval_x = '%sv%s_prop_pval'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         fc_x = '%sv%s_prop_logFC'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         pval_y = '%sv%s_prop_pval'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         fc_y = '%sv%s_prop_logFC'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         st = 'AGG'
#     else:
#         pval_x = '%sv%s_pval'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         fc_x = '%sv%s_logFC'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         pval_y = '%sv%s_pval'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         fc_y = '%sv%s_logFC'%(age_pairs[1][1][0], age_pairs[1][0][0])
#     zscore_x = '%s_zscore'%fc_x
#     zscore_y = '%s_zscore'%fc_y
#     tot = df_table['Tissue'][df_table['Sample.Type']==st].value_counts().reindex(index=tissues)
#     z_pair1_x = (df_table['%s_zscore'%fc_x][df_table['Sample.Type']==st].min(), df_table['%s_zscore'%fc_x][df_table['Sample.Type']==st].max())
#     z_pair2_y = (df_table['%s_zscore'%fc_y][df_table['Sample.Type']==st].min(), df_table['%s_zscore'%fc_y][df_table['Sample.Type']==st].max())
#     if sign == 'up':
#         z_array = np.arange(0, max(z_pair1_x[1], z_pair2_y[1]),0.1)
#     elif sign == 'down':
#         z_array = np.arange(min(z_pair1_x[0], z_pair2_y[0]),0,0.1)
#     else:
#         z_array = np.arange(0, max(abs(z_pair1_x[0]), abs(z_pair1_x[1]), abs(z_pair2_y[0]), abs(z_pair2_y[1])),0.1)
#     x_matrix = np.zeros((len(z_array), len(tissues)+1))
#     y_matrix = np.zeros((len(z_array), len(tissues)+1))
#     x_matrix[:,0] = z_array
#     y_matrix[:,0] = z_array
#     for i, z in enumerate(z_array):
#         if sign == 'up':
#             df_xi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x]>=z)]
#             df_yi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y]>=z)]
#         elif sign == 'down':
#             df_xi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x]<=z)]
#             df_yi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y]<=z)]
#         else:
#             df_xi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x].abs()>=z)]
#             df_yi = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y].abs()>=z)]
#         x_matrix[i,1:] = df_xi['Tissue'].value_counts().reindex(index=tissues)/tot
#         y_matrix[i,1:] = df_yi['Tissue'].value_counts().reindex(index=tissues)/tot
#     df_x = pd.DataFrame(data = x_matrix, columns = ['z']+tissues)
#     df_y = pd.DataFrame(data = y_matrix, columns = ['z']+tissues)
#     df_x.to_csv(os.path.join(fpath_out, '%s_%sv%s_z%s.csv'%(metric, age_pairs[0][1][0], age_pairs[0][0][0], sign)), index = False)
#     df_y.to_csv(os.path.join(fpath_out, '%s_%sv%s_z%s.csv'%(metric,age_pairs[1][1][0], age_pairs[1][0][0], sign)), index = False)
#### part 2: rank age-associated changes comparison between two age pairs, updated on 2020.07.13
# for metric in ['TL', 'AGG', 'Prop']:
#     st = metric
#     if metric == 'Prop':
#         pval_x = '%sv%s_prop_pval'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         fc_x = '%sv%s_prop_logFC'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         pval_y = '%sv%s_prop_pval'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         fc_y = '%sv%s_prop_logFC'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         st = 'AGG'
#     else:
#         pval_x = '%sv%s_pval'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         fc_x = '%sv%s_logFC'%(age_pairs[0][1][0], age_pairs[0][0][0])
#         pval_y = '%sv%s_pval'%(age_pairs[1][1][0], age_pairs[1][0][0])
#         fc_y = '%sv%s_logFC'%(age_pairs[1][1][0], age_pairs[1][0][0])
#     zscore_x = '%s_zscore'%fc_x
#     zscore_y = '%s_zscore'%fc_y
#     if sign == 'up':
#         df_x = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x]>=0)]
#         df_y = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y]>=0)]
#     elif sign == 'down':
#         df_x = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x]<=0)]
#         df_y = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y]<=0)]
#     else:
#         df_x = df_table[(df_table['Sample.Type']==st) & (df_table[pval_x]<=pcutoff) & (df_table[zscore_x].abs()>=0)]
#         df_y = df_table[(df_table['Sample.Type']==st) & (df_table[pval_y]<=pcutoff) & (df_table[zscore_y].abs()>=0)]
#     df_x_out = pd.pivot_table(df_x, values = zscore_x, index = ['N. furzeri Protein Id'], columns = ['Tissue'], aggfunc = np.abs)
#     df_y_out = pd.pivot_table(df_y, values = zscore_y, index = ['N. furzeri Protein Id'], columns = ['Tissue'], aggfunc = np.abs)
#     df_x_out[:] = np.sort(df_x_out.values, axis=0)
#     df_y_out[:] = np.sort(df_y_out.values, axis=0)
#     df_x_out = df_x_out.add_suffix('_%s'%metric)
#     df_y_out = df_y_out.add_suffix('_%s'%metric)
#     df_x_out.to_csv(os.path.join(fpath_out, '%s_%sv%s_absz%s.csv'%(metric, age_pairs[0][1][0], age_pairs[0][0][0], sign)), index = False)
#     df_y_out.to_csv(os.path.join(fpath_out, '%s_%sv%s_absz%s.csv'%(metric,age_pairs[1][1][0], age_pairs[1][0][0], sign)), index = False)
############ End of Section 1: Calculate the extent of changes at increasing z-score cutoff ############


# ############ Section 2: Make the polar heatmap ############  
# #for metric in ['WCL', 'AGG', 'Prop']:
# metric = 'Prop';
# fpath_out = 'Tert/ProliferationIndex'
# pcutoff = 0.05
# age_pairs = [('Young','Old'), ('Old','Tert')]
# # z_cutoff = None#None # 0.75
# tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# sign = 'up'#'both'
# age_comp = '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0])
# # set up the dataframe
# df_TL = pd.read_csv(os.path.join(fpath_out, 'TL_%s_z%s.csv'%(age_comp, sign)))
# df_AGG = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_z%s.csv'%(age_comp, sign)))
# df_Prop = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_z%s.csv'%(age_comp, sign)))
# df_Prop.columns = ['z'] + ['%s_Prop'%i for i in df_Prop.columns[1:]]
# df_heatmap = df_TL.merge(df_AGG, on = 'z', suffixes=('_TL', '_AGG'))
# df_heatmap = df_heatmap.merge(df_Prop, on = 'z')
# df_heatmap.dropna(inplace=True)
# # df_heatmap.dropna(inplace=True, subset=df_heatmap.columns[1:], how='all') 

# n = len(tissues)+1
# m = df_heatmap.shape[0]+1 #indicate number of z-score slices
# rad = np.linspace(0, 10, m)
# a = np.linspace(0, 2 * np.pi, n) + np.pi/6.
# r, th = np.meshgrid(rad, a)

# fig, ax = plt.subplots(1, 3, figsize=(21,6), subplot_kw={"projection":"polar"})
# for i, metric in enumerate(['TL', 'AGG', 'Prop']):
#     z = df_heatmap[['%s_%s'%(tissue, metric) for tissue in tissues]].values.T 
#     im = ax[i].pcolormesh(th, r, z, shading = 'flat', cmap = 'Blues')
#     for a_i in a:
#         ax[i].plot([a_i]*len(rad), rad, ls='--', color = 'lightgrey') 
#     ax[i].set_axis_off()
# plt.colorbar(im)
# plt.show()
# fig.savefig(os.path.join(fpath_out, 'PolarHeatmap%s.pdf'%sign))
# ############ End of Section 2 ############  

# # ############ Section 3: Make wind rose style plot (just pick z>0) to illustrate extent of remodeling ############  
# fpath_out = 'Tert/ProliferationIndex'
# pcutoff = 0.05
# age_pairs = [('Young','Old'), ('Old','Tert')]
# # z_cutoff = None#None # 0.75
# tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594']#, '#3288bd']
# sign = 'both'#'both'
# age_comp = '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0])
# # set up the dataframe
# df_TL = pd.read_csv(os.path.join(fpath_out, 'TL_%s_z%s.csv'%(age_comp, sign)))
# df_AGG = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_z%s.csv'%(age_comp, sign)))
# df_Prop = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_z%s.csv'%(age_comp, sign)))
# df_Prop.columns = ['z'] + ['%s_Prop'%i for i in df_Prop.columns[1:]]
# df_heatmap = df_TL.merge(df_AGG, on = 'z', suffixes=('_TL', '_AGG'))
# df_heatmap = df_heatmap.merge(df_Prop, on = 'z')
# # df_heatmap.dropna(inplace=True)
# # df_heatmap.dropna(inplace=True, subset=df_heatmap.columns[1:], how='all') 

# n = len(tissues)
# m = 5 #range of hexagon traces
# rad = np.linspace(0, 0.2, m)
# a = np.linspace(1, n, n)* 2*np.pi/n - np.pi/6
# # r, th = np.meshgrid(rad, a)

# fig, ax = plt.subplots(1, 3, figsize=(21,6), subplot_kw={"projection":"polar"})
# for i, metric in enumerate(['TL', 'AGG', 'Prop']):
#     z = df_heatmap[['%s_%s'%(tissue, metric) for tissue in tissues]].values.T 
#     for rad_i in rad:
#         ax[i].plot(list(a)+[np.pi/6], [rad_i]*(n+1), ls='--', color = 'lightgrey')
#     for a_r, r in zip(a, z[:,0]):
#         ax[i].plot([0,a_r], [0,0.2], ls='--', color = 'lightgrey') 
#     ax[i].bar(a, z[:,0], color = cm_tissues, edgecolor='k')
#     ax[i].set_title(metric)
#     ax[i].set_axis_off()
# # plt.colorbar(im)
# plt.show()
# fig.savefig(os.path.join(fpath_out, 'WindRose%s.pdf'%sign))
# # ############ End of Section 3 ############  

############ Section 4: Make radar plot (just pick z>0) to illustrate extent of remodeling ############  
# fpath_out = 'Tert/ProliferationIndex'
# pcutoff = 0.05
# age_pairs = [('Young','Old'), ('Old','Tert')]
# # z_cutoff = None#None # 0.75
# tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# sign = 'both'#'both'
# age_comp = '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0])
# # set up the dataframe
# df_TL = pd.read_csv(os.path.join(fpath_out, 'TL_%s_z%s.csv'%(age_comp, sign)))
# df_AGG = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_z%s.csv'%(age_comp, sign)))
# df_Prop = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_z%s.csv'%(age_comp, sign)))
# df_Prop.columns = ['z'] + ['%s_Prop'%i for i in df_Prop.columns[1:]]
# df_heatmap = df_TL.merge(df_AGG, on = 'z', suffixes=('_TL', '_AGG'))
# df_heatmap = df_heatmap.merge(df_Prop, on = 'z')
# # df_heatmap.dropna(inplace=True)
# # df_heatmap.dropna(inplace=True, subset=df_heatmap.columns[1:], how='all') 

# n = len(tissues)
# m = 5 #range of hexagon traces
# rad = np.linspace(0, 0.2, m)# a = np.linspace(1, n, n)* 2*np.pi/n - np.pi/6
# r, th = np.meshgrid(rad, a)

# fig, ax = plt.subplots(1, 3, figsize=(21,6), subplot_kw={"projection":"polar"})
# for i, metric in enumerate(['TL', 'AGG', 'Prop']):
#     z = df_heatmap[['%s_%s'%(tissue, metric) for tissue in tissues]].values.T 
#     ax[i].scatter(a, z[:,0], s= 1000)#shading = 'flat', cmap = 'Blues') 
#     for rad_i in rad:
#         ax[i].plot(list(a)+[np.pi/6], [rad_i]*(n+1), ls='--', color = 'lightgrey')
#     for a_r, r in zip(a, z[:,0]):
#         ax[i].plot([0,a_r], [0,rad[-1]], ls='--', color = 'lightgrey')  
#         ax[i].plot([0,a_r], [0,r], lw=10, color='b')
#     ax[i].set_title(metric)
#     ax[i].set_axis_off()
# # plt.colorbar(im)
# plt.show()
# fig.savefig(os.path.join(fpath_out, 'RadarPlot%s.pdf'%sign))
############ End of Section 3 ############ 

# ############ Section 5: Make stacked histogram to show extent of changes ############ 
# fpath_out = 'Tert'
# comp = 'TvO'
# # set up the dataframe
# df_input = pd.read_csv('kf_combined_prop.csv')
# hist_var = 'logFC_zscore'
# sample_types = ['TL', 'AGG', 'Prop']
# sample_types_markers = ['o', 'o']
# # #age_tag_dict = {'Young': ['128_N', '128_C', '129_N'], 'Old': ['126', '127_N', '127_C'], 'TERT': ['129_C', '130_N', '130_C']}
# # age_tag_dict = {'Young': [1, 2, 3], 'Old': [1, 2, 3], 'TERT': [1, 2, 3]}
# # age_groups = ['Young', 'Old', 'TERT']
# tissues = ['Brain','Gut','Heart','Liver','Muscle','Skin']#, 'Testis']
# tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594']#, '#3288bd']
# for sample_type in sample_types:
#     fig, axes = plt.subplots(len(tissues), 1, sharex=True, sharey=True, figsize=(4,5))
#     for i, tissue in enumerate(tissues):
#         if sample_type !='Prop':
#             df_plot = df_input[(df_input['Sample.Type']==sample_type) & (df_input['Tissue']==tissue)]
#         else:
#             df_plot = df_input[(df_input['Sample.Type']=='AGG') & (df_input['Tissue']==tissue)]
#         hist_coln = '%s_%s'%(comp, hist_var) if sample_type != 'Prop' else '%s_prop_%s'%(comp, hist_var)
#         sns.histplot(data = df_plot, x = hist_coln, kde = True, fill = True, line_kws={ "lw": 1},\
#                   color = tissues_colors[i], ax = axes[i])
#         # axes[i].get_yaxis().set_ticks([])
#         axes[i].set_ylabel(tissue)
#         axes[i].text(0.7, 0.7, 'max = %.2f'%df_plot[hist_coln].abs().max(), ha='left', va='top',transform=axes[i].transAxes)
#         if i == len(tissues) - 1:
#             axes[i].set_xlabel('%s %s %s'%(sample_type, comp, hist_var))
#     axes[i].set_xlim((-6,6))
#     fig.subplots_adjust(hspace=0)
#     plt.savefig(os.path.join(fpath_out, '%s_%s_%s.pdf'%(sample_type, comp, hist_var)),dpi=300)
#     plt.close('all')

# ############ Section 6: Make wind rose style plot (just pick z>0) and add color gradient on z-score ############  
fpath_out = 'Tert/ProliferationIndex'
pcutoff = 0.05
age_pairs = [('Young','Old'), ('Old','Tert')]
# z_cutoff = None#None # 0.75
tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594']#, '#3288bd']
sign = 'both'#'both'
age_comp = '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0])
# set up the dataframe
df_TL = pd.read_csv(os.path.join(fpath_out, 'TL_%s_absz%s.csv'%(age_comp, sign)))
df_AGG = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_absz%s.csv'%(age_comp, sign)))
df_Prop = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_absz%s.csv'%(age_comp, sign)))
df_heatmap = pd.concat([df_TL, df_AGG], axis=1)
df_heatmap = pd.concat([df_heatmap, df_Prop], axis=1)
TL_frac = pd.read_csv(os.path.join(fpath_out, 'TL_%s_z%s.csv'%(age_comp, sign)))
AGG_frac = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_z%s.csv'%(age_comp, sign)))
Prop_frac = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_z%s.csv'%(age_comp, sign)))
Prop_frac.columns = ['z'] + ['%s_Prop'%i for i in Prop_frac.columns[1:]]
df_frac = TL_frac.merge(AGG_frac, on = 'z', suffixes=('_TL', '_AGG'))
df_frac = df_frac.merge(Prop_frac, on = 'z')

n = len(tissues)
m = df_heatmap.shape[0]+1 #indicate number of z-score slices
a = np.linspace(1, n, n)* 2*np.pi/n - np.pi/6

fig, ax = plt.subplots(1, 3, figsize=(21,6), subplot_kw={"projection":"polar"})
for i, metric in enumerate(['TL', 'AGG', 'Prop']):
    z = df_heatmap[['%s_%s'%(tissue, metric) for tissue in tissues]].values.T 
    frac = df_frac[['%s_%s'%(tissue, metric) for tissue in tissues]].values.T 
    for rad_i in a:
        ax[i].plot(list(a)+[np.pi/6], [rad_i/30]*(n+1), ls='--', color = 'lightgrey')
    for a_r, r in zip(a, frac[:,0]):
        ax[i].plot([0,a_r], [0,0.2], ls='--', color = 'lightgrey') 
    ax[i].bar(a, frac[:,0], width = 2*np.pi/n, color = 'None', edgecolor='k')
    for j in range(n):
        z_j = z[j,:][~np.isnan(z[j,:])]
        r_j, th_j = np.meshgrid(np.linspace(0, frac[j,0], z_j.shape[0]), np.linspace(a[j]- np.pi/n, a[(j+1)%n]- np.pi/n + int((j+1)/n)*2*np.pi, 100))
        im = ax[i].pcolormesh(th_j, r_j, np.tile(z_j, (100,1)), cmap = 'Blues')
    ax[i].set_title(metric)
    ax[i].set_axis_off()
plt.colorbar(im)
fig.savefig(os.path.join(fpath_out, 'GradientWindRose%s.pdf'%sign))
# ############ End of Section 3 ############   