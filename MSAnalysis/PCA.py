#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 13:25:56 2019

the file is used to do PCA analysis 

Remove columns that contain "Call" data
Transpose the dataframe so that each row is a sample and each column is a protein
Remove protein description header and set the protein accession numbers as the column headers
Split into train/test sets
Scale values to zero mean and unit varaince if needed
PCA analysis

reading on PCA
https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
https://www.kaggle.com/crawford/principle-component-analysis-gene-expression
https://scikit-learn.org/stable/auto_examples/preprocessing/plot_scaling_importance.html#sphx-glr-auto-examples-preprocessing-plot-scaling-importance-py

@author: yiwenchen

make sure to have input file 'MSResults/kf_combined_prop.csv' ready, change foutpath to any desired output path
"""

import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


import os
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42


def pca_plot(df_geneasrow, feature_key, target_key, out_name, out_dir, title = None, scale = False, n_components = 2, \
             plot = False, target_colors = None, target_shapes = None, legend = False, legend_spacing = 3, \
                 figsize = (16,4), n_vector=2, s = 500, feature_max=10):
#    df_geneasrow.fillna(0,inplace=True)
    df_geneasrow.dropna(inplace=True)
    df = df_geneasrow[target_key].T
    # features
    x = df.values
    # standardize the features
    if scale == True: x = StandardScaler().fit_transform(x)
    # note i am not going to standarize feature because the gene expression is the same type of data
    pca = PCA(n_components=n_components)
    PCs = pca.fit_transform(x) #
    principalDf = pd.DataFrame(data = PCs, columns = ['PC %s'%(i+1) for i in range(n_components)])
    principalDf.insert(0,'target',target_key, True) 
    loadings = pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=df_geneasrow[feature_key])
    loading_matrix = pd.DataFrame(pca.components_.T * np.sqrt(pca.explained_variance_ratio_),\
                                      columns=['PC1', 'PC2'], index=df_geneasrow[feature_key])
    loadings.to_csv(os.path.join(out_dir, '%s_pca_loadings.csv'%out_name))
    loading_matrix.to_csv(os.path.join(out_dir, '%s_pca_standardizedloadings.csv'%out_name))
    top_PC1_vectors = loadings.sort_values(by=['PC1','PC2'], ascending=False).reset_index()
    top_PC2_vectors = loadings.sort_values(by=['PC2','PC1'], ascending=False).reset_index()
    top_PC1_stdv = loading_matrix.sort_values(by=['PC1','PC2'], ascending=False).reset_index()
    top_PC2_stdv = loading_matrix.sort_values(by=['PC2','PC1'], ascending=False).reset_index()
    # visualize
    if plot == True:
        fig, ax = plt.subplots(nrows = 1, ncols=4, figsize = figsize)
        if target_colors == None: target_colors = [matplotlib.cm.get_cmap('Set1')[i] for i in np.linspace(0,1,len(target_key))]
        if target_shapes == None: target_shapes = ['o'] * len(target_key)
        # ax[0] is visualization of sample plotted on PC1 and PC2
        ax[0].set_xlabel('PC 1 (%.2f%%)'%(pca.explained_variance_ratio_[0]*100), fontsize = 15)
        ax[0].set_ylabel('PC 2 (%.2f%%)'%(pca.explained_variance_ratio_[1]*100), fontsize = 15)
        # ax[1] is visualization of sample plotted on PC1 and PC2 as well as the top 5 eigenvector contributing to PC1
        ax[1].set_xlabel('PC 1 (%.2f%%)'%(pca.explained_variance_ratio_[0]*100), fontsize = 15)
        ax[1].set_ylabel('PC 2 (%.2f%%)'%(pca.explained_variance_ratio_[1]*100), fontsize = 15)
        # ax[2] is visualization of sample plotted on PC1 and PC2 as well as the top 5 eigenvector contributing to PC2
        ax[2].set_xlabel('PC 1 (%.2f%%)'%(pca.explained_variance_ratio_[0]*100), fontsize = 15)
        ax[2].set_ylabel('PC 2 (%.2f%%)'%(pca.explained_variance_ratio_[1]*100), fontsize = 15)
        # ax[3] is visualization of all the eigenvector coefficients for PC1 and PC2
        ax[3].set_xlabel('PC 1 (%.2f%%)'%(pca.explained_variance_ratio_[0]*100), fontsize = 15)
        ax[3].set_ylabel('PC 2 (%.2f%%)'%(pca.explained_variance_ratio_[1]*100), fontsize = 15)
        if title ==  None: title  = out_name
        ax[0].set_title(title, fontsize = 15)
        ax[1].set_title('%s PC1 loading'%title, fontsize = 15)
        ax[2].set_title('%s PC2 loading'%title, fontsize = 15)
        ax[3].set_title('%s variable correlation'%title, fontsize = 15)
        for i, target in enumerate(target_key):
            color = target_colors[i]
            shape = target_shapes[i]
            indicesToKeep = (principalDf['target'] == target)
            if (i+1)%legend_spacing == 1:
                ax[0].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = target, edgecolors='black')
                ax[1].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = target, edgecolors='black')
                ax[2].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = target, edgecolors='black')
            else:
                ax[0].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = '', edgecolors='black')
                ax[1].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = '', edgecolors='black')
                ax[2].scatter(principalDf.loc[indicesToKeep, 'PC 1'], principalDf.loc[indicesToKeep, 'PC 2'], \
                       c = color, s = s, marker = shape, label = '', edgecolors='black')
        for j in range(n_vector+1):
            ax[1].arrow(0, 0, top_PC1_vectors.iloc[j,1]*s, top_PC1_vectors.iloc[j,2]*s, color='k', alpha = 0.9, linestyle = '-', linewidth = 1, overhang=0, head_width = 2)
            ax[1].text(top_PC1_vectors.iloc[j,1]*s*1.15,   top_PC1_vectors.iloc[j,2]*s*1.15, top_PC1_vectors.iloc[j,0], fontsize=4)
            ax[1].arrow(0, 0, top_PC1_vectors.iloc[-j,1]*s,top_PC1_vectors.iloc[-j,2]*s, color='k', alpha = 0.9, linestyle = '-', linewidth = 1, overhang=0, head_width = 2)
            ax[1].text(top_PC1_vectors.iloc[-j,1]*s*1.15,  top_PC1_vectors.iloc[-j,2]*s*1.15, top_PC1_vectors.iloc[-j,0], fontsize=4)
            ax[2].arrow(0, 0, top_PC2_vectors.iloc[j,1]*s, top_PC2_vectors.iloc[j,2]*s, color='k', alpha = 0.9, linestyle = '-', linewidth = 1, overhang=0, head_width = 2)
            ax[2].text(top_PC2_vectors.iloc[j,1]*s*1.15,   top_PC2_vectors.iloc[j,2]*s*1.15, top_PC2_vectors.iloc[j,0], fontsize=4)
            ax[2].arrow(0, 0, top_PC2_vectors.iloc[-j,1]*s,top_PC2_vectors.iloc[-j,2]*s, color='k', alpha = 0.9, linestyle = '-', linewidth = 1, overhang=0, head_width = 2)
            ax[2].text(top_PC2_vectors.iloc[-j,1]*s*1.15,  top_PC2_vectors.iloc[-j,2]*s*1.15, top_PC2_vectors.iloc[-j,0], fontsize=4)
        # ax[3] is visualization of all eigenvectors contributing to PC1 and PC2
        ax[3].scatter(top_PC1_stdv.iloc[feature_max:-feature_max,1], top_PC1_stdv.iloc[feature_max:-feature_max, 2], c = 'k', alpha = 0.2, s = 10, edgecolors='black')
        ax[3].scatter(top_PC1_stdv.iloc[:feature_max, 1], top_PC1_stdv.iloc[:feature_max, 2], c = 'r', s = 10, edgecolors='black')
        ax[3].scatter(top_PC1_stdv.iloc[-feature_max:,1], top_PC1_stdv.iloc[-feature_max:, 2], c = 'r', s = 10, edgecolors='black')
        ax[3].scatter(top_PC2_stdv.iloc[:feature_max, 1], top_PC2_stdv.iloc[:feature_max, 2], c = 'r', s = 10, edgecolors='black')
        ax[3].scatter(top_PC2_stdv.iloc[-feature_max:,1], top_PC2_stdv.iloc[-feature_max:, 2], c = 'r', s = 10, edgecolors='black')
        # show legend
        ax[0].legend().set_visible(legend)
        ax[1].legend().set_visible(legend)
        ax[2].legend().set_visible(legend)
        # supress the axis ticklabel because they are pretty much meaningless
        ax[0].set_yticklabels([]);  ax[0].set_xticklabels([])
        # ax[1].set_yticklabels([]);  ax[1].set_xticklabels([]) 
        plt.tight_layout()
        fig.savefig(os.path.join(out_dir, '%s_vector.pdf'%out_name))
        plt.close()
    return pca

# read the data
fpath = ''
fname = 'kf_combined_prop.csv'
pval_cutoff = 0.05
df_all = pd.read_csv(os.path.join(fpath, fname))

# make the long table into short form where the label 
tissues_order = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd'] 
#tissues_order = [ 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', ]
#tissues_colors = ['#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594'] 
age_order = ['Young', 'Old', 'TERT']
sample_types = ['TL', 'AGG']
#age_tag_dict = {'Young': ['Young_X128_C', 'Young_X128_N', 'Young_X129_N'],\
#                'Old': ['Old_X126', 'Old_X127_C', 'Old_X127_N'], \
#                'TERT': ['TERT_X129_C', 'TERT_X130_C', 'TERT_X130_N']}
age_tag_dict = {'Young': ['Young-1','Young-2','Young-3'],\
                'Old': ['Old-1','Old-2','Old-3'],\
                'TERT': ['TERT-1', 'TERT-2', 'TERT-3']}
#age_shape_dict = {'Young':'o', 'Old':'^', 'TERT':'s'}
age_shape_dict = {'Young':'s', 'Old':'o', 'TERT':'^'}
age_color_dict = {'Young':'#91d6e4', 'Old':'#488cca', 'TERT':'#7f7f7f'} # these are the hex color

                  
#### section 1: tissue specific PCA (plot all ages), Figure 1 ####
foutpath = 'PCA'
if not os.path.exists(foutpath): os.makedirs(foutpath)
age_tags = []; age_shapes = []; age_colors = []
for age in age_order:
    for tag in age_tag_dict[age]:
        age_tags.append('log2_%s'%tag)
        age_shapes.append(age_shape_dict[age])
        age_colors.append(age_color_dict[age])
prot_id = 'N. furzeri Protein Id'
# read a tissue at a time, do tissue specific PCA
for s_i, sample_type in enumerate(sample_types):
    for t_i, tissue in enumerate(tissues_order):
        age_tags_i = age_tags[:-3] if tissue == 'Testis' else age_tags
        df_tissue = df_all[['Human', prot_id]+age_tags_i][(df_all['Tissue']==tissue) & (df_all['Sample.Type'] == sample_type)]
        pca_i = pca_plot(df_tissue, prot_id, age_tags_i, '%s_%s_pca_nolegend'%(sample_type, tissue), foutpath,\
                  title ='%s %s'%(sample_type, tissue), scale = True, plot = True, target_colors = age_colors, target_shapes = age_shapes)
        new_coln = ['%s_%s'%(tissue, coln) if 'log2' in coln else coln for coln in df_tissue.columns]
        df_tissue.columns = new_coln
        if t_i == 0:
            df_sample = df_tissue.iloc[:,1:].copy(deep=True)
        else:
            df_sample = df_sample.merge(df_tissue.iloc[:,1:], left_on = prot_id, right_on = prot_id, how = 'outer')
    sample_target_key = df_sample.columns[1:]
    pca_s = pca_plot(df_sample, prot_id, sample_target_key, '%s_pca_nolgend'%sample_type, foutpath, title = sample_type, figsize = (16,4),\
              scale = True, plot = True, legend = False, target_colors = [i for i in tissues_colors for j in range(len(age_tags))],\
              target_shapes = age_shapes*len(tissues_order), legend_spacing=9)
#### End of section 1: tissue specific PCA (plot all ages) ####


#### section 2: tissue specific PCA (plot just old and young) ####
foutpath = 'PCA/OYOnly'
if not os.path.exists(foutpath): os.makedirs(foutpath)
age_order = ['Young', 'Old']
age_tags = []; age_shapes = []; age_colors = []
for age in age_order:
    for tag in age_tag_dict[age]:
        age_tags.append('log2_%s'%tag)
        age_shapes.append(age_shape_dict[age])
        age_colors.append(age_color_dict[age])
prot_id = 'N. furzeri Protein Id'
# read a tissue at a time, do tissue specific PCA
for s_i, sample_type in enumerate(sample_types):
    for t_i, tissue in enumerate(tissues_order):
        df_tissue = df_all[['Human', prot_id]+age_tags][(df_all['Tissue']==tissue) & (df_all['Sample.Type'] == sample_type)]
        pca_i = pca_plot(df_tissue, prot_id, age_tags, '%s_%s_pca_nolegend'%(sample_type, tissue), foutpath, figsize = (16,4),\
                title ='%s %s'%(sample_type, tissue), scale = True, plot = True, target_colors = age_colors, target_shapes = age_shapes)
        new_coln = ['%s_%s'%(tissue, coln) if 'log2' in coln else coln for coln in df_tissue.columns]
        df_tissue.columns = new_coln
        # uncomment this section if you desire to plot all tissue in one plot
        if t_i == 0:
            df_sample = df_tissue.iloc[:,1:].copy(deep=True)
        else:
            df_sample = df_sample.merge(df_tissue.iloc[:,1:], left_on = prot_id, right_on = prot_id, how = 'outer')
    sample_target_key = df_sample.columns[1:]
    pca_s = pca_plot(df_sample, prot_id, sample_target_key, '%s_pca_nolgend'%sample_type, foutpath, title = sample_type, figsize = (16,4),\
            scale = True, plot = True, legend = False, target_colors = [i for i in tissues_colors for j in range(len(age_tags))],\
            target_shapes = age_shapes*len(tissues_order), legend_spacing=9)
#### End of section 2: tissue specific PCA (only old and young data) ####
