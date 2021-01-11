#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 14:39:13 2019
Updated on 2020.05.13

@author: yiwenchen
"""

#import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib
import seaborn as sns
import matplotlib.backends.backend_pdf


def tag_age(tag):
    if 'Young' in tag:
        out = 'Young'
    elif 'Old' in tag:
        out = 'Old'
    elif 'TERT' in tag:
        out = 'TERT'
    return out

# function used to plot barplot with jitter for categorial data (per tissue)
# function used to plot barplot with jitter for categorial data (per tissue)
def bar_sample(df, protein_name, identifier, st, tissues=[], fname_out='', fpath_out='',\
               ax='', ind = True, ages = ['Young','Old'], ctype='default', ylabel=True):#, 'Tert']):
#    age_tag_dict = {'Young': ['Young_X128_C', 'Young_X128_N', 'Young_X129_N'],\
#                'Old': ['Old_X126', 'Old_X127_C', 'Old_X127_N'], \
#                'Tert': ['TERT_X129_C', 'TERT_X130_C', 'TERT_X130_N']}
    age_tag_dict = {'Young': ['Young-1','Young-2','Young-3'],\
                'Old': ['Old-1','Old-2','Old-3'],\
                'TERT': ['TERT-1', 'TERT-2', 'TERT-3']}
    age_colors = {'Young':'#91d6e4', 'Old':'#488cca', 'TERT':'#7f7f7f'} # ['#99d8c9','#e5f5f9']
    age_alphas = {'Young':0.5, 'Old':1, 'TERT':1}
    tissue_colors = {'Brain':'#d53e4f', 'Gut':'#fc8d59', 'Heart':'#fee08b',\
                     'Liver':'#ffffbf', 'Muscle':'#e6f598', 'Skin':'#99d594', 'Testis':'#3288bd'}
    if st != 'Prop':
        value_tag = ['log2_%s'%i for age in ages for i in age_tag_dict[age]]
        value_coln = '%s abundance'%st
        st_fil = st
    else:
        value_tag = ['log2_%s_prop'%i for age in ages for i in age_tag_dict[age]]
        value_coln = 'aggregation propensity'
        st_fil = 'AGG'
    pid = df['N. furzeri Final Symbol'][df[identifier]==protein_name].unique()
    if len(pid)== 0: 
        print ("protein isn't in df")
        return None, None
    pid = pid[0]
    if len(tissues)==0: tissues = df['Tissue'][df[identifier]==protein_name].unique()
    df_data = pd.melt(df[['Tissue']+value_tag][(df[identifier]==protein_name)&(df['Tissue'].isin(tissues))&(df['Sample.Type']==st_fil)],\
                      id_vars=['Tissue'],var_name='age',value_name=value_coln)
    df_data['age'] = df_data['age'].apply(tag_age)
    if (len(tissues)>1) & (len(tissues)>len(df_data['Tissue'].unique())):
        for i in (set(tissues) - set(df_data['Tissue'].unique())):
            df_add = pd.DataFrame({'Tissue':[i]*len(ages)*3, 'age':ages*3})
            df_data = df_data.append(df_add, ignore_index=True)
#    df_data.sort_values(by=['Tissue','age'], ascending=[True,False], inplace=True)
    df_data.sort_values(by=['Tissue'], ascending=[True], inplace=True)
    if df_data.empty == True: 
        print ("check to make sure you specify the tissues where the protein is expressed in")
        return None, None
    if ax == '': fig, ax = plt.subplots(figsize = (6,8)) if len(tissues) == 1 else plt.subplots(figsize = (len(tissues)*4,10))
    if ax != '': fig = None
    ######## plot with catplot ##########
    sns.axes_style("white")
    if len(tissues) == 1:
        # plot the aggregate portion
        if ctype == 'default':
            sns.boxplot(ax = ax, x = 'age', y = value_coln, data = df_data, \
            palette = {'Young': tissue_colors[tissues[0]],'Old':tissue_colors[tissues[0]], 'TERT':age_colors['TERT']}, linewidth=2.5) 
        else:
            sns.boxplot(ax = ax, x = 'age', y = value_coln, data = df_data, palette = age_colors, linewidth=2.5)
        sns.stripplot(ax = ax, x = 'age', y = value_coln, data = df_data, jitter = True, dodge = True, size = 24, color = 'black')
        if ctype == 'default':
            for patch, age in zip(ax.artists, ages):
                r, g, b, a = patch.get_facecolor()
                patch.set_facecolor((r, g, b, age_alphas[age]))     
        if ylabel == False: ax.yaxis.label.set_visible(False)                                                           
        ax.set_title('%s: %s'%(tissues[0], pid), fontsize=12)
    else:
        # plot the aggregate portion
        sns.boxplot(ax=ax, x='Tissue', y=value_coln, hue='age', data=df_data, hue_order = ages, palette=age_colors, linewidth=2.5)           #sns.color_palette()          
        sns.stripplot(ax=ax, x='Tissue', y=value_coln, hue='age', data=df_data, hue_order = ages, jitter=True, dodge=True, size=24, color='black')                                            
        ax.set_title(pid,fontsize=12)
    ax.xaxis.label.set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    if len(ages) == 3:
        ax.legend(handles[0:3], labels[0:3], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    elif len(ages) == 2:
        ax.legend(handles[0:2], labels[0:2], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    ######################## End of option 1 #####################
    # save the figure
    if ind == True:
        if fname_out == '':
            fname_out = '%s_%s_%s.pdf'%(st, pid.upper(), tissues[0]) if len(tissues)==1 else '%s_%s.pdf'%(st, pid.upper())
        plt.savefig(os.path.join(fpath_out, fname_out), bbox_inches = 'tight')
#        plt.tight_layout()
    else:
        if fig != None:
            if fname_out.startswith('XP_'):
                fig.suptitle(fname_out[:-4])
            else:
                fig.suptitle(fname_out)
    plt.close()
    return fig, ax


def bar_MS(df, protein_name, identifier, tissues=[], fname_out='', fpath_out='',ind = True, ages = ['Young','Old']):#, 'TERT']):
    pid = df['N. furzeri Final Symbol'][df[identifier]==protein_name].unique()
    if len(pid)== 0:  return None, None
    pid = pid[0]
    if len(tissues) == 0: tissues = df['Tissue'][df[identifier]==protein_name].unique()
    if len(tissues) == 1:
        fig, axes = plt.subplots(1,2, figsize = (10,10)) 
    else:
        fig, axes = plt.subplots(2,1, figsize =(len(tissues)*6,20))    
    fig_tl,  ax_tl = bar_sample(df, protein_name, identifier, 'TL', tissues=tissues, ind = False, ax = axes.flat[0])
    fig_agg, ax_agg = bar_sample(df, protein_name, identifier, 'AGG', tissues=tissues, ind = False, ax = axes.flat[1])
    if ind == True:
        if fname_out == '': fname_out = 'TLAGG_%s.pdf'%pid.upper()
        fig.savefig(os.path.join(fpath_out, fname_out), bbox_inches = 'tight')
    else:
        if fname_out.startswith('XP_'):
            fig.suptitle(fname_out[:-4])
        else:
            fig.suptitle(fname_out)
    plt.close()
    return fig, axes


def bar_MS_separate(df, protein_name, identifier, tissues=[], fname_out='', fpath_out='',\
                    ind = True, ages = ['Old','TERT'], stype = ['TL','AGG', 'Prop'], highlight = [], sharey=True):#, 'TERT']):
    pid = df['N. furzeri Final Symbol'][df[identifier]==protein_name].unique()
    if len(pid)== 0: return None, None
    pid = pid[0]
    if len(tissues) == 0: tissues = df['Tissue'][df[identifier]==protein_name].unique()
    fig, axes = plt.subplots(len(stype),len(tissues), figsize =(len(tissues)*6,len(stype)*8), sharey=sharey)
    for  i, st in enumerate(stype):
        for j, tissue in enumerate(tissues): 
            ctype='default' if tissue in highlight else 'other'
            ylabel = True if j == 0 else False
            if len(stype) == 1:
                ax = axes[j]
            elif len(tissues) == 1:
                ax = axes[i]
            else:
                ax = axes[i,j]
            ax.set_ylabel(st)
            fig_ij, ax_ij = bar_sample(df, protein_name, identifier, st, tissues=[tissue],\
                                       ages = ages, ind = False, ax = ax, ylabel = ylabel, ctype= ctype)
    if ind == True:
        if fname_out == '': fname_out = '%s_%s.pdf'%(''.join(stype), pid.upper())
        fig.savefig(os.path.join(fpath_out, fname_out), bbox_inches = 'tight')
    else:
        if fname_out.startswith('XP_'):
            fig.suptitle(fname_out[:-4])
        else:
            fig.suptitle(fname_out, fontsize=10)
    plt.close()
    return fig, axes


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


def combined_scatter(df, age1, age2, prot_id, id_key, fpath_out = '', tissues='all', ax=None, lim=None):
    tissue_colors_dict = {'Brain':'#d53e4f', 'Gut':'#fc8d59', 'Heart':'#fee08b',\
                          'Liver':'#ffffbf', 'Muscle':'#e6f598', 'Skin':'#99d594', 'Testis':'#3288bd'}
#    age_tag_dict = {'Young': ['Young_X128_C', 'Young_X128_N', 'Young_X129_N'],\
#            'Old': ['Old_X126', 'Old_X127_C', 'Old_X127_N'], \
#            'TERT': ['TERT_X129_C', 'TERT_X130_C', 'TERT_X130_N']}
    age_tag_dict = {'Young': ['Young-1','Young-2','Young-3'],\
                'Old': ['Old-1','Old-2','Old-3'],\
                'TERT': ['TERT-1', 'TERT-2', 'TERT-3']}
    out_coln = ['log2_%s'%i for age in [age1,age2] for i in age_tag_dict[age]]
    df_agg = df[['Tissue']+out_coln][(df[id_key]==prot_id) & (df['Sample.Type']=='AGG')]
    df_tl = df[['Tissue']+out_coln][(df[id_key]==prot_id) & (df['Sample.Type']=='TL')]
    human = df['Human'][df[id_key]==prot_id].unique()[0]
    pid = df['N. furzeri Protein Id'][df[id_key]==prot_id].unique()[0]
    if (df_tl.empty != True) & (df_agg.empty != True):
        df_agg = pd.melt(df_agg,id_vars=['Tissue'],var_name='Age', value_name='AGG')
        df_tl = pd.melt(df_tl,id_vars=['Tissue'],var_name='Age', value_name='TL')
        df_plot = df_agg.merge(df_tl, on=['Tissue','Age'], how = 'outer')
        df_plot['Age'] = df_plot['Age'].apply(lambda x: age1 if age1 in x else age2)
        if ax == None:
            fig, ax = plt.subplots(figsize = (5,4))
            sns.set(font_scale = 1) 
            sns.set_style("white")
            sns.scatterplot(x='TL', y='AGG', data=df_plot, hue='Age', style='Tissue', ax=ax)
            ax.set_title(human)
            ax.legend(markerscale=1, bbox_to_anchor=(1,1))
            plt.tight_layout()
            fig.savefig(os.path.join(fpath_out, '%s_%s_scatter.pdf'%(human, pid)))
            plt.close(fig)
    return ax


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 24})
#sns.set(font_scale = 5) # use 5 for paper figure
sns.set_style("white")


#### Section 1: Given a list, output panel of old and young AGG or TL ####
df_all = pd.read_csv(os.path.join('','kf_combined_prop.csv'))
fig_path = 'ProteinBoxPlot'
if not os.path.exists(fig_path): os.makedirs(fig_path)
protein = 'cct4'; tissue = 'Gut'; identifier = 'N. furzeri Final Symbol'
# ##### part 1: plot only TL or AGG
bar_sample(df_all, protein, identifier, 'AGG', tissues=[tissue], fpath_out=fig_path, ages = ['Young','Old', 'TERT'])
# ##### part 2: plot all MS data with shared y-axis
bar_MS(df_all, protein, identifier, fpath_out=fig_path)
# ##### part 3: plot data from selective tissues with where each tissue is a separate subplot
tissues=['Gut', 'Liver']
fig, ax = bar_MS_separate(df_all, protein, identifier, stype = ['TL','AGG','Prop'], tissues=tissues,\
                highlight=tissues, fpath_out=fig_path, sharey = False)
########################### End of Section 2 ###########################
