#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:18:24 2019

# this is the script used to create the raw table with functional information on brain prop hits

@author: yiwenchen
"""

import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats
import itertools
from collections import OrderedDict 
#import xlsxwriter


mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 2)
sns.set_style("white", {'axes.linewidth':0.5})


# define a function to get stats on properties of interest when meeting certain cutoff
def filter_stats(df, tissues, var, cutoff, direction = '>='):
    out = []
    pos_out = []
    tot_out = []
    for tissue in tissues:
        df_fil = df[['N. furzeri Protein Id',var]][df['Tissue']==tissue].dropna(axis=0)
        if direction == '>=':
            pos = df_fil['N. furzeri Protein Id'][(df_fil[var]>=cutoff)].count()
        elif direction == '<=':
            pos = df_fil['N. furzeri Protein Id'][(df_fil[var]<=cutoff)].count()
        elif direction == '>':
            pos = df_fil['N. furzeri Protein Id'][(df_fil[var]>cutoff)].count()
        elif direction == '<':
            pos = df_fil['N. furzeri Protein Id'][(df_fil[var]<cutoff)].count()
        tot = df_fil['N. furzeri Protein Id'].count()
        out.append(pos*1./tot)
        tot_out.append(tot)
        pos_out.append(pos)
        print ('in %s, pos=%s, tot = %s, frac = %s'%(tissue, pos, tot, pos*1./tot))
    return out, pos_out, tot_out


# define a function to get stats on properties of interest when meeting certain cutoff
def filter_list(df, tissues, df_list, id_coln = 'Associated genes'):
    out = []
    pos_out = []
    tot_out = []
    for tissue in tissues:
        df_fil = df[['N. furzeri Protein Id','Human']][df['Tissue']==tissue].dropna(axis=0)
        pos = df_fil['N. furzeri Protein Id'][df_fil['Human'].isin(df_list[id_coln].tolist())].count()
#        print df_fil['N. furzeri Protein Id'][df_fil['Human'].isin(df_list[id_coln].tolist())]
        tot = df_fil['N. furzeri Protein Id'].count()
        out.append(pos*1./tot)
        tot_out.append(tot)
        pos_out.append(pos)
        print ('in %s, pos=%s, tot = %s, frac = %s'%(tissue, pos, tot, pos*1./tot))
    return out, pos_out, tot_out


def fisher_exact(row):
    oddsratio, pvalue = scipy.stats.fisher_exact([[row['Cutoff Y All'], row['Cutoff Y Hit']], [row['Cutoff N All'], row['Cutoff N Hit']]])
    return pvalue


def age_comp(df, age1, age2, st, pval_cutoff = 0.05, fc = 'both'):
    if st == 'Prop':
#        fc_coln = 'Propensity_%s-%s_FC'%(age1, age2)
#        pval_coln = 'Propensity_%s-%s_pval'%(age1, age2)
        fc_coln = '%sv%s_prop_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
        df_data = df[df['Sample.Type']=='AGG'].sort_values(['Tissue', fc_coln], ascending = [True, False])
    else:
        fc_coln = '%sv%s_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_pval'%(age1[0], age2[0])
        df_data = df[df['Sample.Type']==st].sort_values(['Tissue', fc_coln], ascending = [True, False]) 
    out_coln = ['Tissue','Sample.Type','N. furzeri symbol','N. furzeri Protein Id','N. furzeri Final Symbol','vertebrate only', '# Vertebrates with orthologs (total 39)','Human',fc_coln, pval_coln]
    out_coln = [i for i in out_coln if i in df.columns]
    if fc == 'both':
        df_out = df_data[out_coln][df_data[pval_coln]<=pval_cutoff]
        df_out.loc[:,'Hit.Type'] = '%s_%s'%(st, pval_cutoff)
    elif fc == 'up':
        df_out = df_data[out_coln][(df_data[pval_coln]<=pval_cutoff) & (df_data[fc_coln]>1)]
        df_out.loc[:,'Hit.Type'] = '%s_%s_%s'%(st, pval_cutoff, fc) 
    elif fc == 'down':
        df_out = df_data[out_coln][(df_data[pval_coln]<=pval_cutoff) & (df_data[fc_coln]<1)]
        df_out.loc[:,'Hit.Type'] = '%s_%s_%s'%(st, pval_cutoff, fc)        
    return df_out


def sig_fracc(df, age1, age2, st, pval_cutoff = 0.05, fc_cutoff = 1, direction = 'both', table = True): # fc_cutoff = 0.2 means age2_vs_age1 has 20% increase or decrease
    if st == 'Prop':
#        fc_coln = 'Propensity_%s-%s_FC'%(age1, age2)
#        pval_coln = 'Propensity_%s-%s_pval'%(age1, age2)
        fc_coln = '%sv%s_prop_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
        df_data = df[df['Sample.Type']=='AGG'].sort_values(['Tissue', fc_coln], ascending = [True, False])
    else:
        fc_coln = '%sv%s_FC'%(age1[0], age2[0])
        pval_coln = '%sv%s_pval'%(age1[0], age2[0])
        df_data = df[df['Sample.Type']==st].sort_values(['Tissue', fc_coln], ascending = [True, False]) 
    out_coln = ['Tissue','Sample.Type','N. furzeri symbol','N. furzeri Protein Id','N. furzeri Final Symbol','vertebrate only', '# Vertebrates with orthologs (total 39)','Human',fc_coln, pval_coln]
    if direction == 'both':
        df_out = df_data[(df_data[fc_coln]<=(1.-fc_cutoff))|(df_data[fc_coln]>=(1+fc_cutoff))]
        df_out = df_out[out_coln][df_out[pval_coln]<=pval_cutoff]
        if df_out.empty == False: df_out.loc[:,'Hit.Type'] = '%s_%s'%(st, pval_cutoff)
    elif direction == 'up':
        df_out = df_data[out_coln][(df_data[pval_coln]<=pval_cutoff) & (df_data[fc_coln]>=(1.+fc_cutoff))]
        if df_out.empty == False: df_out.loc[:,'Hit.Type'] = '%s_%s_%s'%(st, pval_cutoff, direction) 
    elif direction == 'down':
        df_out = df_data[out_coln][(df_data[pval_coln]<=pval_cutoff) & (df_data[fc_coln]<=(1.-fc_cutoff))]
        if df_out.empty == False: df_out.loc[:,'Hit.Type'] = '%s_%s_%s'%(st, pval_cutoff, direction)        
    if df_out.empty == False:
        df_out_stats = pd.pivot_table(df_out, values=fc_coln, index=['Tissue', 'Sample.Type'], aggfunc='count')
    else:
        df_out_stats = None 
    return df_out, df_out_stats


################## Section 1: Make clean table that contain hit of interest and their type info ##########
# look through matrix then loop through different age comparisons
# either with or without z_cutoff #
### part 1: age-associated changes, updated on 2020.04.20        
hit_path = 'Correlomics/Comp'
if not os.path.exists(hit_path): os.makedirs(hit_path)
pcutoff = 0.05
z_cutoff = None#None # 0.75
df_table = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
for metric in ['TL', 'AGG', 'Prop']: 
    st = metric
    for age1, age2 in itertools.combinations(['Young', 'Old', 'TERT'], 2):
        out_name = '%sAll_%sv%s'%(metric, age2[0], age1[0])
        out_coln = ['Tissue', 'N. furzeri Final Symbol', 'N. furzeri Protein Id', 'Human', 'log2_%s'%age1, 'log2_%s'%age2]
        if metric == 'Prop':
            pval = '%sv%s_prop_pval'%(age2[0], age1[0])
            fc = '%sv%s_prop_logFC'%(age2[0], age1[0])
            st = 'AGG'
            out_coln += [fc, '%sv%s_logFC'%(age2[0], age1[0]), 'log2_Young_TL']
        else:
            pval = '%sv%s_pval'%(age2[0], age1[0])
            fc = '%sv%s_logFC'%(age2[0], age1[0])
            out_coln += [fc, '%sv%s_prop_logFC'%(age2[0], age1[0]), 'log2_Young_TL']
        df_out_all = df_table[out_coln][df_table['Sample.Type']==st].copy()
        df_out_all.loc[:,'Hit Type'] = out_name
        df_out_all.to_csv(os.path.join(hit_path, '%s.csv'%out_name), index = False)
        if z_cutoff == None:
            df_out_sig = df_table[out_coln][(df_table['Sample.Type']==st) & (df_table[pval]<=pcutoff)].copy()
            df_out_sig.loc[:,'Hit Type'] = np.where(df_out_sig[fc]>=0, out_name.replace('All', 'SigPos'), out_name.replace('All', 'SigNeg'))
            df_out_sig.to_csv(os.path.join(hit_path, '%s.csv'%(out_name.replace('All', 'Sig'))), index = False)
        else:
            zscore = '%s_zscore'%fc
            df_out_sig = df_table[out_coln][(df_table['Sample.Type']==st) & (df_table[pval]<=pcutoff) & (df_table[zscore].abs()>=z_cutoff)].copy()
            df_out_sig.loc[:,'Hit Type'] = np.where(df_out_sig[fc]>=0, out_name.replace('All', 'SigPosZ%s'%z_cutoff), out_name.replace('All', 'SigNegZ%s'%z_cutoff))
            df_out_sig.to_csv(os.path.join(hit_path, '%s.csv'%(out_name.replace('All', 'SigZ%s'%z_cutoff))), index = False)
### part 2: protein with high AGG or Prop for specific age, updated on 2020.04.30 
hit_path = 'Correlomics/AgeSpecific'
if not os.path.exists(hit_path): os.makedirs(hit_path)
z_cutoff = 2
direction = 'greater'
df_table = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
for metric in ['TL', 'AGG', 'Prop']: 
    for age in ['Young', 'Old', 'TERT']:
        out_name = '%s_%sZ%s'%(age, metric, z_cutoff)
        out_coln = ['Tissue', 'N. furzeri Final Symbol', 'N. furzeri Protein Id', 'Human']
        if metric == 'Prop':
            value = 'log2_%s_prop'%age
            df_out_all = df_table[df_table['Sample.Type']=='AGG'].copy()
        else:
            value = 'log2_%s'%age
            df_out_all = df_table[df_table['Sample.Type']==metric].copy()
        df_out_all.loc[:,'Hit Type'] = out_name
        zscore = '%s_zscore'%value
        out_coln += [value, zscore]
        df_out_sig = df_out_all[out_coln][df_out_all[zscore].abs()>=z_cutoff] if direction == 'greater' else df_out_all[out_coln][df_out_all[zscore].abs()<=z_cutoff].copy()
        df_out_sig.loc[:,'Hit Type'] = np.where(df_out_sig[zscore]>=0, '%sPosZ%s'%(metric, z_cutoff), '%sNegZ%s'%(metric, z_cutoff))
        df_out_sig.to_csv(os.path.join(hit_path, '%s.csv'%out_name), index = False)            
################################# End of Section 1 ##############################


################# Section 3: Output TL, AGG, AND Prop Hits #################
#go through matrix then loop through different age comparisons
in_path = ''
out_path = 'SigChanges'
if not os.path.exists(out_path): os.makedirs(out_path)
fname = 'kf_combined_prop.csv'
df_all = pd.read_csv(os.path.join(in_path, fname))
sts = ['TL', 'AGG', 'Prop']
age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
pval_cutoff = 0.05
for age_pair in age_pairs:
    age1, age2 = age_pair
    for st in sts:
        df_out = age_comp(df_all, age1, age2, st, pval_cutoff = pval_cutoff)
        df_out.to_csv(os.path.join(out_path, '%s_%sv%s.csv'%(st, age1[0], age2[0])), index = False)
############################### End of Section 3 ##############################


################ Section 4: Look at PrD enrichment in hits and make bar plots #################
in_path = ''
in_fname = 'kf_combined_prop.csv'
map_path = 'Properties'
map_fname = 'Killifish_DISOLLRCider.csv'
out_path = 'SigChanges'
if not os.path.exists(out_path): os.makedirs(out_path)
df_all = pd.read_csv(os.path.join(in_path, in_fname))
df_map = pd.read_csv(os.path.join(map_path, map_fname))
df_all = df_all.merge(df_map, left_on =['N. furzeri Gene Id','N. furzeri Protein Id'] , right_on = ['Gene Id','N. furzeri Protein Id'], how='left')
sts = ['AGG']
var = 'LLR'
var_cutoff = 0
zcutoff = 0
age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
tissues = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis'] 
bar_palette = ['#dc9fdc', '#D3D3D3']#9300d2']
for st in sts:
    for age_pair in age_pairs:
        age1, age2 = age_pair
        if st == 'Prop':
            fc_coln = '%sv%s_prop_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']=='AGG'].sort_values(['Tissue', fc_coln], ascending = [True, False])
        else:
            fc_coln = '%sv%s_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']==st].sort_values(['Tissue', fc_coln], ascending = [True, False])
        z_coln = '%s_zscore'%fc_coln
        all_out, all_pos, all_tot = filter_stats(df_data, tissues, var, var_cutoff, direction = '>')
        hit_out, hit_pos, hit_tot = filter_stats(df_data[(df_data[fc_coln]>=0) & (df_data[pval_coln]<=0.05) & (df_data[z_coln].abs()>=zcutoff)], tissues, var, var_cutoff, direction = '>')
        data_bar = {'Tissue':tissues*2,'Type':['PosSig']*len(tissues)+['Detected']*len(tissues),\
                    'Frac with PrD': hit_out + all_out, 'with PrD': hit_pos + all_pos,\
                        'tot':hit_tot + all_tot}
        data_test = {'Tissue':tissues,'Cutoff Y All':all_pos, 'Cutoff Y Hit':hit_pos,\
                      'Cutoff N All':np.asarray(all_tot)-np.asarray(all_pos), 'Cutoff N Hit':np.asarray(hit_tot)-np.asarray(hit_pos)}
        df_bar = pd.DataFrame.from_dict(data_bar, orient='columns')
        df_test = pd.DataFrame.from_dict(data_test)
        df_test['Fisher pval'] = df_test.apply(fisher_exact, axis=1)
        if zcutoff == 0:
            df_test.to_csv(os.path.join(out_path, '%s_%sv%s_Pos%s%s.csv'%(st, age1[0], age2[0],var, var_cutoff)))
        else:
            df_test.to_csv(os.path.join(out_path, '%s_%sv%s_Pos%s%s_z%s.csv'%(st, age1[0], age2[0],var, var_cutoff, zcutoff)))
        fig, ax = plt.subplots(figsize=(10,6))
        sns.barplot(x='Tissue', y='Frac with PrD', hue='Type', data=df_bar, ax = ax, palette = bar_palette, edgecolor='.1')
        ax.set_title('%s %sv%s %s>%s'%(st, age1[0], age2[0], var, var_cutoff))
        if zcutoff == 0:
            fig.savefig(os.path.join(out_path, '%s_%sv%s_Pos%s%s.pdf'%(st, age1[0], age2[0],var, var_cutoff)))
        else:
            fig.savefig(os.path.join(out_path, '%s_%sv%s_Pos%s%s_z%s.pdf'%(st, age1[0], age2[0],var, var_cutoff, zcutoff)))
############################## End of Section 4 ##############################


#################### Section 5: Look at IDR enrichment in hits and make bar plots #################        
in_path = ''
in_fname = 'kf_combined_prop.csv'
map_path = 'Properties'
map_fname = 'Killifish_DISOLLRCider.csv'
out_path = 'SigChanges'
if not os.path.exists(out_path): os.makedirs(out_path)
df_all = pd.read_csv(os.path.join(in_path, in_fname))
df_map = pd.read_csv(os.path.join(map_path, map_fname))
df_all = df_all.merge(df_map, left_on =['N. furzeri Gene Id','N. furzeri Protein Id'] , right_on = ['Gene Id','N. furzeri Protein Id'], how='left')
sts = ['AGG']
var = 'Disordered Frac'
var_cutoff = 0.3
zcutoff = 0
age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
tissues = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis'] 
bar_palette = ['#dc9fdc', '#D3D3D3']#9300d2']
for st in sts:
    for age_pair in age_pairs:
        age1, age2 = age_pair
        if st == 'Prop':
            fc_coln = '%sv%s_prop_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']=='AGG'].sort_values(['Tissue', fc_coln], ascending = [True, False])
        else:
            fc_coln = '%sv%s_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']==st].sort_values(['Tissue', fc_coln], ascending = [True, False])
        z_coln = '%s_zscore'%fc_coln
        all_out, all_pos, all_tot = filter_stats(df_data, tissues, var, var_cutoff, direction = '>')
        hit_out, hit_pos, hit_tot = filter_stats(df_data[(df_data[fc_coln]>=0) & (df_data[pval_coln]<=0.05) & (df_data[z_coln].abs()>=zcutoff)], tissues, var, var_cutoff, direction = '>=')
        data_bar = {'Tissue':tissues*2, 'Type':['PosSig']*len(tissues)+['Detected']*len(tissues), 'Frac with IDR': hit_out + all_out}
        data_test = {'Tissue':tissues, 'Cutoff Y All':all_pos, 'Cutoff Y Hit':hit_pos,\
                      'Cutoff N All':np.asarray(all_tot)-np.asarray(all_pos), 'Cutoff N Hit':np.asarray(hit_tot)-np.asarray(hit_pos)}
        df_bar = pd.DataFrame.from_dict(data_bar)
        df_test = pd.DataFrame.from_dict(data_test)
        df_test['Fisher pval'] = df_test.apply(fisher_exact, axis=1)
        if zcutoff == 0:
            df_test.to_csv(os.path.join(out_path,'%s_%sv%s_Diso%s.csv'%(st, age1[0], age2[0], var_cutoff)))
        else:
            df_test.to_csv(os.path.join(out_path,'%s_%sv%s_Diso%s_z%s.csv'%(st, age1[0], age2[0], var_cutoff, zcutoff)))
        fig, ax = plt.subplots(figsize=(10,6))
        sns.barplot(x='Tissue', y='Frac with IDR', hue='Type', data=df_bar, ax = ax, palette = bar_palette, edgecolor='.1')
        ax.set_title('%s %sv%s'%(st, age1[0], age2[0]))
        if zcutoff == 0:
            fig.savefig(os.path.join(out_path, '%s_%sv%s_Diso%s.pdf'%(st, age1[0], age2[0], var_cutoff)))
        else:
            fig.savefig(os.path.join(out_path, '%s_%sv%s_Diso%s_z%s.pdf'%(st, age1[0], age2[0], var_cutoff, zcutoff)))
################################## End of Section 5 ##############################
        
        
##################### Section 6: Look at ND enrichment in hits and make bar plots #################        
in_path = ''
in_fname = 'kf_combined_prop.csv'
map_path = 'BrainFig'
map_fname = '20201115_CuratedND.csv'
out_path = 'BrainFig/NDEnrichment'
if not os.path.exists(out_path): os.makedirs(out_path)
df_all = pd.read_csv(os.path.join(in_path, in_fname))
df_map = pd.read_csv(os.path.join(map_path, map_fname))
sts = ['TL']
zcutoff = 0
age_pairs = [('Old', 'Young'), ( 'Tert', 'Young'), ('Tert', 'Old')]
tissues = [ 'Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis'] 
bar_palette = ['#dc9fdc', '#D3D3D3']#9300d2']
for st in sts:
    for age_pair in age_pairs:
        age1, age2 = age_pair
        if st == 'Prop':
            fc_coln = '%sv%s_prop_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_prop_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']=='AGG'].sort_values(['Tissue', fc_coln], ascending = [True, False])
        else:
            fc_coln = '%sv%s_logFC'%(age1[0], age2[0])
            pval_coln = '%sv%s_pval'%(age1[0], age2[0])
            df_data = df_all[df_all['Sample.Type']==st].sort_values(['Tissue', fc_coln], ascending = [True, False])
        z_coln = '%s_zscore'%fc_coln
        all_out, all_pos, all_tot = filter_list(df_data, tissues, df_map)
        hit_out, hit_pos, hit_tot = filter_list(df_data[(df_data[pval_coln]<=0.05) & (df_data[z_coln]>=zcutoff)], tissues, df_map)
        data_bar = {'Tissue':tissues*2,'Type':['PosSig']*len(tissues)+['Detected']*len(tissues), 'ND Association':hit_out + all_out}
        data_test = {'Tissue':tissues,'Cutoff Y All':all_pos, 'Cutoff Y Hit':hit_pos, \
                      'Cutoff N All': np.asarray(all_tot)-np.asarray(all_pos), 'Cutoff N Hit':np.asarray(hit_tot)-np.asarray(hit_pos)}
        df_bar = pd.DataFrame.from_dict(data_bar)
        df_test = pd.DataFrame.from_dict(data_test)
        df_test['Fisher pval'] = df_test.apply(fisher_exact, axis=1)
        zsuffix = '' if zcutoff == 0 else '_z%s'%zcutoff
        df_test.to_csv(os.path.join(out_path,'%s_%sv%s_ND%s.csv'%(st, age1[0], age2[0], zsuffix)), index=False)
        fig, ax = plt.subplots(figsize=(10,6))
        sns.barplot(x='Tissue', y='ND Association', hue='Type', data=df_bar, ax = ax, palette = bar_palette, edgecolor='.1')
        ax.set_title('%s %sv%s'%(st, age1[0], age2[0]))
        fig.savefig(os.path.join(out_path, '%s_%sv%s_ND%s.pdf'%(st, age1[0], age2[0], zsuffix)))
################################### End of Section 6 ##############################
