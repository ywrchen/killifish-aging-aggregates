#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 20:21:04 2019

this is the script to generate reproducibility plots in a super compact format for supp figure

@author: yiwenchen
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import scipy.stats
import seaborn as sns
import itertools 

#configure the default font and style
mpl.rcParams['pdf.fonttype'] = 42
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Helvetica']

sns.set(style="white")

# the simplest figure generator where the pair of values are provided to generate panels of figures
def mbyn_subplots(df, pair_text, pair_labels, row_num, col_num, title, logscale, na_method, xlim, ylim, diagonal_line, directory):
    if len(pair_labels) > (row_num*col_num):
        sys.exit('make sure you allocate enough space for all subplots')
    fig, axes = plt.subplots(row_num, col_num, sharex='col', sharey='row', figsize=(6,6))
    if row_num!=1:
        ax_list = [item for sublist in axes for item in sublist] 
    else:
        ax_list = axes
    for ax, labels, ax_text in zip(ax_list, pair_labels, pair_text):
        column_x, column_y = labels
        if na_method == 'drop':
            df_nonull = df[[column_x, column_y]].dropna()
        elif na_method == 'fill':
            df_nonull = df[[column_x, column_y]].fillna()
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
            ax.annotate('pearsonr = %.2f'%(coeff), xy = (0.05, 0.9), xycoords='axes fraction')
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
    fig.text(0.5, 0.04, 'normalized peptide spectra counts', ha='center')
    fig.text(0.04, 0.5, 'normalized peptide spectra counts', va='center', rotation='vertical')
    plt.suptitle(title)
    plt.savefig(os.path.join(directory, '%s.pdf'%title.replace(" ","_")))
    plt.close('all')
    return 0

# Visualize variations in biological replicates (ratio, PSMs, PSMsNorm, agg/wcl_PSMsNorm) and output the all figures in a big pannel where each row represent each data from one tissue

def replicates(df, tissue, sample_type, age_groups, replicate_num, variable, logscale=1, na_method = 'drop', xlim = [10e-3, 10e3], ylim = [10e-3, 10e3], diagonal_line = True, directory = ''):
    pair_labels = []
    pair_text = []
    row_num = len(age_groups)
    replicate_num
    for age in age_groups:
        tags = ['%s-%s'%(age,i) for i in range(1, replicate_num+1)]
        labels = ['%s_%s: %s_%s'%(sample_type, tissue, i, variable) for i in tags]
        pair_labels+= [[labels[0], labels[1]], [labels[0], labels[2]], [labels[1], labels[2]]]
        pair_text+= ['%s vs. %s'%(tags[1], tags[0]), '%s vs. %s'%(tags[2], tags[0]), '%s vs. %s'%(tags[2], tags[1])]
    title = '%s %s %s'%(tissue, sample_type, variable)
    mbyn_subplots(df, pair_text, pair_labels, row_num, replicate_num, title, logscale, na_method, xlim, ylim, diagonal_line, directory)
    return 0


############ Section 1: 7 by 1 grid , Figure 1  ############
# Generate a panel of plots where the logInt of one biological replicate is plotted against another
sample_types = ['AGG']
#age_tag_dict = {'Young': ['128_N', '128_C', '129_N'], 'Old': ['126', '127_N', '127_C'], 'TERT': ['129_C', '130_N', '130_C']}
age_tag_dict = {'Young': [1, 2, 3], 'Old': [1, 2, 3], 'TERT': [1, 2, 3]}
age_groups = ['Young']
tissues = ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']

input_fname = 'kf_combined_prop.csv'
input_fpath = ''
out_fpath = 'Reproducibility'
if not os.path.exists(out_fpath): os.makedirs(out_fpath)
df_input = pd.read_csv(os.path.join(input_fpath, input_fname))

row_comb = [('%s-%s'%(age, comb[0]),'%s-%s'%(age, comb[1])) for age in age_groups for comb in itertools.combinations(age_tag_dict[age],2)]
row_num = len(row_comb)*len(sample_types)
col_num = len(tissues)
fig, axes = plt.subplots(row_num, col_num, sharex=True, sharey=True, figsize=(25,20))
for i, row_i in enumerate(row_comb):
    for i_add, sample_type in enumerate(sample_types):
        for j, tissue in enumerate(tissues):
            if (tissue == 'Testis') & ('TERT' in row_i[0]):
                continue
            else:
                ax_iter = axes[i*len(sample_types)+i_add,j]
                df_iter = df_input[(df_input['Sample.Type']==sample_type) & (df_input['Tissue']==tissue)]
                coln_x = 'log2_%s'%row_i[0]
                coln_y = 'log2_%s'%row_i[1]
                # calculate the linear regression
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df_iter[coln_x], df_iter[coln_y])
                # calculate the pearson correlation coefficients
                coeff, p_value = scipy.stats.pearsonr(df_iter[coln_x], df_iter[coln_y])
                # plot scatter plot
                ax_iter.scatter(df_iter[coln_x], df_iter[coln_y], c='k', alpha = 0.2)
                ax_iter.set_title('%s %s'%(tissue, sample_type), size = 'small')
                ax_iter.annotate('pearson r = %.2f'%(coeff), xy = (0.05, 0.9), xycoords='axes fraction')
                # plot x = y line
                lims = [np.min([ax_iter.get_xlim(), ax_iter.get_ylim()]), np.max([ax_iter.get_xlim(), ax_iter.get_ylim()])]
                # now plot both limits against eachother
                ax_iter.plot(lims, lims, color = '#d8dcd6', ls='--')
                ax_iter.set_aspect('equal')
                ax_iter.set_xlabel(coln_x)
                ax_iter.set_ylabel(coln_y)
                ax_iter.set_xlim(lims)
                ax_iter.set_ylim(lims)
mpl.rcParams.update({'font.size': 20})
plt.tight_layout()
plt.savefig(os.path.join(out_fpath, 'Reproducibility_YoungAGG.pdf'))
plt.close('all')
############ Section 1 ############


#################################
############ Section 2: simplied 14 by 9 grid with limited axis label, Figure S1 #############
# Generate a panel of plots where the logInt of one biological replicate is plotted against another    
sample_types = ['TL', 'AGG']
sample_types_markers = ['o', 'o']
#age_tag_dict = {'Young': ['128_N', '128_C', '129_N'], 'Old': ['126', '127_N', '127_C'], 'TERT': ['129_C', '130_N', '130_C']}
age_tag_dict = {'Young': [1, 2, 3], 'Old': [1, 2, 3], 'TERT': [1, 2, 3]}
age_groups = ['Young', 'Old', 'TERT']
tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin', 'Testis']
tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd'] 

input_fname = 'kf_combined_prop.csv'
input_fpath = ''
out_fpath = 'Reproducibility'
if not os.path.exists(out_fpath): os.makedirs(out_fpath)
df_input = pd.read_csv(os.path.join(input_fpath, input_fname))

lims = (-2,10)
#row_comb = [('%s_X%s'%(age, comb[0]),'%s_X%s'%(age, comb[1])) for age in age_groups for comb in itertools.combinations(age_tag_dict[age],2)]
row_comb = [('%s-%s'%(age, comb[0]),'%s-%s'%(age, comb[1])) for age in age_groups for comb in itertools.combinations(age_tag_dict[age],2)]
text_comb = [(rep_no[0],'%s_%s'%(rep_no[1], age)) for age in age_groups for rep_no in itertools.combinations([1, 2, 3],2)]
row_num = len(row_comb)
col_num = len(tissues)*len(sample_types) 
fig, axes = plt.subplots(row_num, col_num, sharex=True, sharey=True, figsize=(28,18))
stats_dict = dict()
for i, row_i in enumerate(row_comb):
    for j_add, sample_type in enumerate(sample_types):
        for j, tissue in enumerate(tissues):
            if (tissue == 'Testis') & ('TERT' in row_i[0]):
                axes[i,j*2+j_add].set_axis_off()
            else:
                ax_iter = axes[i,j*2+j_add]
                df_iter = df_input[(df_input['Sample.Type']==sample_type) & (df_input['Tissue']==tissue)]
                coln_x = 'log2_%s'%row_i[0]
                coln_y = 'log2_%s'%row_i[1]
                # calculate the linear regression
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df_iter[coln_x], df_iter[coln_y])
                # calculate the pearson correlation coefficients
                coeff, p_value = scipy.stats.pearsonr(df_iter[coln_x], df_iter[coln_y])
                # plot scatter plot
                ax_iter.scatter(df_iter[coln_x], df_iter[coln_y], c=tissues_colors[j], marker = sample_types_markers[j_add], alpha = 0.2)
                age = text_comb[i][1][2:]; rep_x = text_comb[i][0]; rep_y = text_comb[i][1][0]
                text1 = 'TL' if sample_type == 'WCL' else sample_type
                if i == 0:
                    ax_iter.set_title('%s %s'%(text1, age), fontsize = 22)
                if j*2+ j_add ==0:
                    ax_iter.set_ylabel('%s vs %s '%(rep_y, rep_x), fontsize = 24)
                ax_iter.annotate('r = %.2f'%(coeff), xy = (0.2, 0.05), xycoords='axes fraction', fontsize = 24)
                stats_dict['(%s, %s)'%(i,(j*2+j_add))] = [tissue, sample_type, age, 'rep%svs%s'%(rep_y, rep_x), slope, intercept, coeff, p_value]
                # now plot both limits against eachother
                ax_iter.plot(lims, lims, color = '#d8dcd6', ls='--')
                ax_iter.set_aspect('equal')
                if j*2 + j_add == 0: 
                    la = 1
                ax_iter.set_xlim(lims)
                ax_iter.set_ylim(lims)
                ax_iter.tick_params(axis='x', labelsize=16)
                ax_iter.tick_params(axis='y', labelsize=16)
fig.subplots_adjust(wspace=0, hspace=0)
plt.savefig(os.path.join(out_fpath, 'AllReproducibility_long_multicolor.png'),dpi=300)
plt.close('all')  
# output the correlation to csv
df_stats = pd.DataFrame.from_dict(stats_dict, orient='index', columns = ['tissue', 'sample.type','age',\
           'yvsx', 'linear.regression.slope','linear.regression.intercept','pearson.r','pearson.pval'])      
df_stats.to_csv(os.path.join(out_fpath, 'AllReproducibility.csv'))
############ End of Section 2: 4 by 9 grid #### ########
