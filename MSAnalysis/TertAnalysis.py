#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 23:07:52 2019

this is used to filter out and generate venn diagram for protein
that are aggregating both in old an tert compared to young

@author: yiwenchen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats
import matplotlib
import seaborn as sns
import matplotlib.gridspec as gridspec

matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams.update({'font.size': 22})


def filter_data(df, st, filter_col, filter_val, pval_coln, sign = 'pos'):
    tissues = set(df['Tissue'].tolist())
    tissue_ind = dict()
    for tissue in tissues:
        if sign == 'pos':
            tissue_ind[tissue] = df['N. furzeri Protein Id'][(df['Sample.Type'] == st) & (df['Tissue'] == tissue) & (df[pval_coln] <= 0.05) & (df[filter_col] >= filter_val)].tolist()
        else:
            tissue_ind[tissue] = df['N. furzeri Protein Id'][(df['Sample.Type'] == st) & (df['Tissue'] == tissue) & (df[pval_coln] <= 0.05) & (df[filter_col] <= filter_val)].tolist()
    return tissue_ind


def r(x, y):
    return scipy.stats.pearsonr(x, y)[0]


class SeabornFig2Grid():

    def __init__(self, seaborngrid, fig,  subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
            isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n,m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i,j], self.subgrid[i,j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h= self.sg.ax_joint.get_position().height
        h2= self.sg.ax_marg_x.get_position().height
        r = int(np.round(h/h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r+1,r+1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        #https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure=self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())


def col_names(st, age_pair):
    if st == 'Prop':
        st_out='AGG'
        fc_coln = 'Delta_Propensity_%s-%s'%(age_pair[0], age_pair[1])
        pval_coln = 'Propensity_%s-%s_pval'%(age_pair[0], age_pair[1])
    else:
        st_out = st
        fc_coln = '%sv%s_logFC'%(age_pair[0][0], age_pair[1][0])
        pval_coln = '%sv%s_pval'%(age_pair[0][0], age_pair[1][0])
    return fc_coln, pval_coln, st_out


def age_color(age):
    age_colors = ['#91d6e4','#488cca', '#7f7f7f']
    if 'Young' in age: 
        c = age_colors[0]
    elif 'Old' in age: 
        c = age_colors[1]
    elif 'Tert' in age:
        c = age_colors[2]
    return c


# ########### Section 13: Compare extend of change across tissues for OvY and TvO ###########
# # look through matrix then loop through different age comparisons
# # either with or without z_cutoff #
# fpath_out = 'Tert/ProliferationIndex'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# pcutoff = 0.05
# age_pairs = [('Young','Old'), ('Old','Tert')]
# z_cutoff = None#None # 0.75
# df_table = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
# tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# sign = 'both'
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
# ## part 2: visualize the results from part1 updated on 2020.06.01 (show both age-pairs)
# metric = 'Prop';
# age_comp = ('%sv%s'%(age_pairs[0][1][0], age_pairs[0][0][0]), '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0]))
# df_pair1 = pd.read_csv(os.path.join(fpath_out, '%s_%s_z%s.csv'%(metric, age_comp[0], sign)))
# df_pair2 = pd.read_csv(os.path.join(fpath_out, '%s_%s_z%s.csv'%(metric, age_comp[1], sign)))
# df_heatmap = df_pair1.merge(df_pair2, on = 'z', suffixes=('_%s'%age_comp[0], '_%s'%age_comp[1]))
# new_order = ['z']+['%s_%s'%(i,j) for i in tissues for j in age_comp]
# df_heatmap =  df_heatmap[new_order]
# for i in range(len(tissues)):
#     df_heatmap.insert(i*3+1, 'blank_%s'%i, np.nan)
# f, ax = plt.subplots(figsize=(8,6))
# # Draw the heatmap with the mask and correct aspect ratio
# sns.set(font_scale = 1)
# matrix_map = sns.heatmap(df_heatmap.iloc[:,2:], cmap='Blues',vmin=0, vmax=0.17,\
#                         square=False, ax = ax, cbar_kws={"shrink": .5}, yticklabels = False)
# matrix_map.set_facecolor('Gray')
# matrix_map.get_figure().savefig(os.path.join(fpath_out, '%s_%sv%s_z%s.pdf'%(metric, age_comp[0], age_comp[1], sign)))
# # part 3: visualize the results from part1 updated on 2020.06.02, show just TvO results but combine the tissue types
# age_comp = '%sv%s'%(age_pairs[1][1][0], age_pairs[1][0][0])
# df_TL = pd.read_csv(os.path.join(fpath_out, 'TL_%s_z%s.csv'%(age_comp, sign)))
# df_AGG = pd.read_csv(os.path.join(fpath_out, 'AGG_%s_z%s.csv'%(age_comp, sign)))
# df_Prop = pd.read_csv(os.path.join(fpath_out, 'Prop_%s_z%s.csv'%(age_comp, sign)))
# df_Prop.columns = ['z'] + ['%s_Prop'%i for i in df_Prop.columns[1:]]
# df_heatmap = df_TL.merge(df_AGG, on = 'z', suffixes=('_TL', '_AGG'))
# df_heatmap = df_heatmap.merge(df_Prop, on = 'z')
# df_heatmap.dropna(inplace=True)
# # df_heatmap.dropna(inplace=True,subset=df_heatmap.columns[1:], how='all')
# df_heatmap.insert(len(tissues)+1, 'blank_1', np.nan)
# df_heatmap.insert(len(tissues)*2+2, 'blank_2', np.nan)
# f, ax = plt.subplots(figsize=(8,6))
# # Draw the heatmap with the mask and correct aspect ratio
# sns.set(font_scale = 1)
# matrix_map = sns.heatmap(df_heatmap.iloc[:,1:], cmap='Blues',vmin=0, vmax=0.17,\
#                         square=False, ax = ax, cbar_kws={"shrink": .5}, yticklabels = False)
# matrix_map.set_facecolor('Gray')
# matrix_map.get_figure().savefig(os.path.join(fpath_out, 'All_%s_z%s.pdf'%(age_comp[1], sign)))
# ################################## End of Section 13 ##############################


########### Section 14: Correlate cell cycle status to extend of change across tissues (TvO) ###########
fpath_data = 'Tert/ProliferationIndex'
df_all = pd.read_csv(os.path.join(fpath_data, 'FUCCI_YWT.csv'))
df_data = df_all.iloc[:-1,:]
# #### section 1: correlation between extend of AGG or PROP and proliferation index
fig, ax = plt.subplots(1,3, figsize=(12, 4))
x0 = 'late mitosis'; y0 = 'G1'; s0 = 'extend of AGG'; s1 = 'extend of PROP'
color_dict = {'Brain':'#d53e4f', 'Gut':'#fc8d59', 'Heart':'#fee08b',\
                          'Liver':'#ffffbf', 'Muscle':'#e6f598', 'Skin':'#99d594'}#, 'Testis':'#3288bd'}
slope_0, intercept_0, r_value_0, p_value_0, std_err_0 = scipy.stats.linregress(df_data[x0],df_data[y0])
sns.regplot(x=x0, y=y0, data=df_data, ax=ax[0], ci=None, color= 'grey', \
            scatter_kws={'s': df_data[s0]*1000, 'color':list(color_dict.values()), 'edgecolors':'k'},\
            line_kws={'label':"$r^2$={0:.2f}".format(r_value_0**2), 'color': 'grey'});
ax[0].set_xlim(0,100)
ax[0].legend()
slope_10, intercept_10, r_value_10, p_value_10, std_err_10 = scipy.stats.linregress(df_data[s0],df_data[x0])
slope_11, intercept_11, r_value_11, p_value_11, std_err_11 = scipy.stats.linregress(df_data[s0],df_data[y0])
# sns.regplot(x=s0, y=x0, data=df_data, ax=ax[1], color= 'grey',\
#             scatter_kws={'s': df_data[y0]*10, 'color':list(color_dict.values()), 'edgecolors':'k'},\
#             ci=None, line_kws={'label': x0 +" $r^2$={0:.2f}".format(r_value_10**2), 'color': 'grey', 'alpha':0.5});
sns.regplot(x=s0, y=y0, data=df_data, ax=ax[1], color='red', \
            scatter_kws={'s': df_data[x0], 'color':list(color_dict.values()), 'edgecolors':'k'},\
            ci=None, line_kws={'label': y0 +" $r^2$={0:.2f}".format(r_value_11**2), 'color': 'red', 'alpha':0.5});
ax[1].set_xlim(0,0.18)
ax[1].set_ylim(0,50)
ax[1].legend()
slope_20, intercept_20, r_value_20, p_value_20, std_err_20 = scipy.stats.linregress(df_data[s1],df_data[x0])
slope_21, intercept_21, r_value_21, p_value_21, std_err_21 = scipy.stats.linregress(df_data[s1],df_data[y0])
# sns.regplot(x=s1, y=x0, data=df_data, ax=ax[2], color= 'grey',\
#             scatter_kws={'s': df_data[y0]*10, 'color':list(color_dict.values()), 'edgecolors':'k'},\
#             ci=None, line_kws={'label': x0 +" $r^2$={0:.2f}".format(r_value_20**2), 'color': 'grey', 'alpha':0.5});
sns.regplot(x=s1, y=y0, data=df_data, ax=ax[2], color='red', \
            scatter_kws={'s': df_data[x0], 'color':list(color_dict.values()), 'edgecolors':'k'},\
            ci=None, line_kws={'label': y0 +" $r^2$={0:.2f}".format(r_value_21**2), 'color': 'red', 'alpha':0.5});
ax[2].set_xlim(0,0.15)
ax[2].set_ylim(0,50)
ax[2].legend()
plt.tight_layout()
fig.savefig(os.path.join(fpath_data, 'FUCCI_Y_G1.pdf'))
#### section 2: make stacked bar plot ####
fig2, ax2 = plt.subplots(1,1, figsize=(4, 4))
df_data.iloc[:,:-2].plot(x = 'tissue', kind = 'barh', stacked = True, mark_right = True, ax = ax2,\
                        color = ['black','red','green','orange'])
fig2.savefig(os.path.join(fpath_data,'FUCCI_YWT_stackedbar.pdf')) 