#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 17:00:00 2019

@author: yiwenchen
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy.stats
import matplotlib
from os import listdir
from os.path import isfile, join
import seaborn as sns

##configure the default font and style
#from matplotlib import rcParams
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Helvetica']
##plt.style.use('ggplot')
#
#rcParams.update({'font.size': 9})
##plt.rcParams['svg.fonttype'] = 'none'

# the simplest figure generator where the pair of values are provided to generate panels of figures
def mbyn_subplots(df, pair_text, pair_labels, row_num, col_num, title, logscale, na_method, xlim, ylim, diagonal_line, directory, axis_label):
    if len(pair_labels) > (row_num*col_num):
        sys.exit('make sure you allocate enough space for all subplots')
    fig, axes = plt.subplots(row_num, col_num, sharex='col', sharey='row', figsize=(6,6))
    if row_num!=1:
        ax_list = [item for sublist in axes for item in sublist] 
    else:
        ax_list = axes
    for i, df_i, labels, ax_text in zip(range(len(pair_labels)), df, pair_labels, pair_text):
        ax = ax_list[i]
        column_x, column_y = labels
        if na_method == 'drop':
            df_nonull = df_i[[column_x, column_y]].dropna()
        elif na_method == 'fill':
            df_nonull = df_i[[column_x, column_y]].fillna()
        matrix_row, matrix_coln = df_nonull.shape
        #print df_nonull.shape
        if matrix_row <=1:
            print ('%s and %s has only one entry, so no scatter_hist_joint plot will be generated.' %(column_x, column_y))
        else:
            # calculate the linear regression
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df_nonull[column_x], df_nonull[column_y])
            # calculate the pearson correlation coefficients
            coeff, p_value = scipy.stats.pearsonr(df_nonull[column_x], df_nonull[column_y])
            # calculate the pearson correlation coefficients on the log-transformed data
            coeff_ln, p_value_ln = scipy.stats.pearsonr(np.log(df_nonull[column_x]), np.log(df_nonull[column_y]))
            # plot scatter plot
            ax.scatter(df_nonull[column_x], df_nonull[column_y], c='k', alpha = 0.2)
            ax.set_title(ax_text, size = 'small')
            #ax.annotate('pearsonr = %.2f, ln_r=%.2f'%(coeff,coeff_ln), xy = (0.05, 0.9), xycoords='axes fraction')
            #ax.annotate('pr=%.2f,ln_r=%.2f'%(coeff,coeff_ln), xy = (0.05, 0.9), xycoords='axes fraction')
            ax.annotate('pearsonr = %.2f'%(coeff_ln), xy = (0.05, 0.9), xycoords='axes fraction')
            if logscale == 1:
                ax.set_xscale('log')
                ax.set_yscale('log')
            if logscale == 2:
                ax.set_xscale('log',basex = 2)
                ax.set_yscale('log',basey = 2)
            if xlim != [np.nan, np.nan]:
                ax.set_xlim(xlim[0],xlim[1])
            if ylim != [np.nan, np.nan]:
                ax.set_ylim(ylim[0],ylim[1])
            if diagonal_line == True:
                # plot x = y line
                lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
                # now plot both limits against eachother
                ax.plot(lims, lims, color = '#d8dcd6', ls='--', c='0.3')
                #    ax.set_aspect('equal')
                ax.set_xlim(lims)
                ax.set_ylim(lims)
    while i+1 < len(ax_list):
        i+=1
        ax_list[i].axis('off')
    fig.text(0.5, 0.04, axis_label[0], ha='center')
    fig.text(0.04, 0.5, axis_label[1], va='center', rotation='vertical')
    plt.suptitle(title)
    plt.savefig(os.path.join(directory, '%s.pdf'%title.replace(" ","_").replace('/','')))
    plt.close('all')
    return 0


def new_rank(df, fc = 'logFC', pval = 'Qval'):
    df.loc[:,'rankstats'] = df[fc]*(-1)*np.log10(df[pval])
    return df


############## Section 1: run this if you want to create sample_type specific summary table on protein of interest ##############
# #### updated on 2020.04.27
# # geneset_category is the main vary to specificy which gene set to analyze
# geneset_cateogry = 'Ribosome' ##Chaperones','Proteasome','Ribosome','UbiquitinProteolysis','Autophagy', 'Lysosome'
# # sample_type can be TL, or AGG, or Prop
# sample_type = 'TL' #AGG, TL, Prop
# pairs = [('Young','Old'),('Young','Tert'),('Old','Tert')]
# fpath = 'GeneSets/Killifish/%s'%geneset_cateogry
# foutpath = 'GeneSets/MS/%s_%s'%(geneset_cateogry, sample_type)
# if not os.path.exists(foutpath): os.makedirs(foutpath)
# geneset_fnames = [f for f in listdir(fpath) if isfile(join(fpath, f))]
# geneset_fnames = [i for i in geneset_fnames if '_Killifish.csv' in i]
# all_geneset = []
# geneset_list = []
# df_input = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
# #geneset_fnames = ['Other_chaperones_Killifish.csv']
# for geneset_fname in geneset_fnames:
#     with open(os.path.join(fpath, geneset_fname), 'r') as f: geneset = f.read()
#     geneset = geneset.replace('\r','\n')
#     geneset = geneset.replace('\n\n','\n')
#     geneset = list(set(geneset.split('\n')))
#     all_geneset += geneset
#     geneset_list.append(geneset)
# for pair in pairs:
#     age_pair = '%sv%s'%(pair[1][0], pair[0][0])
#     for geneset, set_name, i in zip(geneset_list, geneset_fnames,range(len(geneset_list))):
#         #### generate table that's convenient for heatmap or clustering ####
#         if sample_type == 'Prop':
#             df_data = df_input[(df_input['N. furzeri Final Symbol'].isin(geneset)) & (df_input['Sample.Type']=='AGG')]
#             fc = '%s_prop_logFC'%age_pair
#             pval = '%s_prop_pval'%age_pair
#         else:
#             df_data = df_input[(df_input['N. furzeri Final Symbol'].isin(geneset)) & (df_input['Sample.Type']==sample_type)]
#             fc = '%s_logFC'%age_pair
#             pval = '%s_pval'%age_pair
#         if (df_data[fc].dropna(how='all').empty) == False:
#             df_data = new_rank(df_data, fc = fc, pval = pval)
#             df_out = df_data.pivot_table(index = 'Human', columns = 'Tissue', values = [fc, pval, 'rankstats'])
#             df_out.insert(0, 'classification', set_name.replace('_Killifish.csv','').replace('_All',''))
#             df_out_flat = df_data.loc[:,['Human', 'Tissue', fc, pval, 'rankstats']]
#             df_out_flat.insert(0, 'classification', set_name.replace('_Killifish.csv','').replace('_All',''))
#             df_out_flat.to_csv(os.path.join(foutpath, set_name.replace('Killifish.csv', '%s_%s_flat.csv'%(age_pair,sample_type))),index = False)
#             if i == 0:
#                 df_allout = df_out.copy(deep = True)
#                 df_allout_flat = df_out_flat.copy(deep = True)
#             else:
#                 df_allout = df_allout.append(df_out, ignore_index = False)
#                 df_allout_flat = df_allout_flat.append(df_out_flat, ignore_index = False)
#             out_name = set_name.replace('_Killifish', '_%s_%s'%(age_pair,sample_type))
#             df_out.to_csv(os.path.join(foutpath, out_name),index = True)
#     if len(geneset_fnames) >1:
#         df_allout.to_csv(os.path.join(foutpath, 'All%s_%s_%s.csv'%(fpath[fpath.rindex('/')+1:], age_pair, sample_type)),index = True)
#         df_allout_flat.to_csv(os.path.join(foutpath, 'All%s_%s_%s_flat.csv'%(fpath[fpath.rindex('/')+1:], age_pair, sample_type)),index = False)
############# End of Section 1 ##############


############### Section 2: run this if you want to create sample_type specific summary table on protein of interest ##############
####### the output is used for the Cytoscape visualization
#### updated on 2020.04.27
# geneset_cateogries = ['Autophagy','Chaperones','Lysosome','Proteasome','Ribosome','UbiquitinProteolysis']#'Proteasome' ##Chaperones'
# sample_type = 'AGG'
# pair = ('Young','Old')
# age_pair = '%sv%s'%(pair[1][0], pair[0][0])
# #pairs = [('Young','Old'),('Young','Tert'),('Old','Tert')]
# df_input = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
# for set_i, geneset_cateogry in enumerate(geneset_cateogries):
#     fpath = 'GeneSets/Killifish/%s'%geneset_cateogry
#     foutpath = 'GeneSets/MS'#%(geneset_cateogry, sample_type)
#     if not os.path.exists(foutpath): os.makedirs(foutpath)
#     geneset_fnames = [f for f in listdir(fpath) if isfile(join(fpath, f))]
#     geneset_fnames = [i for i in geneset_fnames if '_Killifish.csv' in i]
#     geneset_list = []
#     #geneset_fnames = ['Other_chaperones_Killifish.csv']
#     for geneset_fname in geneset_fnames:
#         with open(os.path.join(fpath, geneset_fname), 'r') as f: geneset = f.read()
#         geneset = geneset.replace('\r','\n')
#         geneset = geneset.replace('\n\n','\n')
#         geneset = list(set(geneset.split('\n')))
# #        all_geneset += geneset
#         geneset_list.append(geneset)
#     for geneset, set_name, i in zip(geneset_list, geneset_fnames, range(len(geneset_list))):
#         #### generate table that's convenient for heatmap or clustering ####
#         if sample_type == 'Prop':
#             df_data = df_input[(df_input['N. furzeri Final Symbol'].isin(geneset)) & (df_input['Sample.Type']=='AGG')]
#             fc = '%s_prop_logFC'%age_pair
#             pval = '%s_prop_pval'%age_pair
#         else:
#             df_data = df_input[(df_input['N. furzeri Final Symbol'].isin(geneset)) & (df_input['Sample.Type']==sample_type)]
#             fc = '%s_logFC'%age_pair
#             pval = '%s_pval'%age_pair
#         if (df_data[fc].dropna(how='all').empty) == False:
#             df_data = new_rank(df_data, fc = fc, pval = pval)
#             df_out_flat = df_data.loc[:,['Human', 'Tissue', fc, pval, 'rankstats']]
#             df_out_flat.insert(0, 'classification', set_name.replace('_Killifish.csv','').replace('_All',''))
#             if set_i == 0:
#                 df_allout_flat = df_out_flat.copy(deep = True)
#             else:
#                 df_allout_flat = df_allout_flat.append(df_out_flat, ignore_index = False)
# df_allout_flat.to_csv(os.path.join(foutpath, 'AllGenesets_%s_%s_flat.csv'%(age_pair, sample_type)),index = False)
############### End of Section 2 ##############
