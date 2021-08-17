#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 18:04:10 2021

purpose: scan killifish proteome and look protein with KFERQ-containing peptides
The motif is always flanked by
a glutamine on one of the sides (as the pentapeptide functions as a targeting sequence in both directions) 
and contains one or two of the positive residues K and R, 
one or two of the hydrophobic residues F, L, I or V and 
one of the negatively charged E or D residues

@author: yiwenchen
"""

import os
import pandas as pd
import numpy as np
import re
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import matplotlib as mpl



mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 20})
sns.set(style="white")


#Functions
# Read a fasta file and output a proteome_dict = {proteinA:proteinA_sequence, proteinB:proteinB_sequence...}
# File_Name = '2016.02.12 - nfurzeri_proteins_final_mline_symbol.fasta'
def read_old_fasta(File_Name, motif = None):
    f = open(File_Name)
    Raw_Data = f.read()
    f.close()
    Raw_Data = Raw_Data.split('>')[1:]
    Len_Prot = len(Raw_Data)
    Proteome_Dict = dict()
    for i in np.arange(0, Len_Prot):
        Entry = Raw_Data[i]
        Prot_Last_Index = Entry.find('\n')
        Protein = Entry[:Prot_Last_Index].replace(' ','')
        Sequence = Entry[Prot_Last_Index:].replace('\n','').replace(' ' ,'')
        if motif == None:
            Proteome_Dict[Protein] = str(Sequence)
        else:
            for m_i in motif:
                regexp = re.compile(m_i)
                if regexp.search(Sequence):
                    Proteome_Dict[Protein] = str(Sequence)
    return Proteome_Dict


def new_rank(df, fc = 'logFC', pval = 'Qval'):
    df['rankstats'] = df[fc]*(-1)*np.log10(df[pval])
    return df


# ########## Section #1: Look for KFERQ-like peptide in killifish proteome ##########
# motif = ['Q[(K,R)]{1,2}[(F,L,I,V)]{1,2}[E|D]', '[(K,R)]{1,2}[(F,L,I,V)]{1,2}[E|D]Q',\
#          'Q[(K,R)]{1,2}[E|D][(F,L,I,V)]{1,2}', '[(K,R)]{1,2}[E|D][(F,L,I,V)]{1,2}Q',\
#          'Q[(F,L,I,V)]{1,2}[(K,R)]{1,2}[E|D]', '[(F,L,I,V)]{1,2}[(K,R)]{1,2}[E|D]Q',\
#          'Q[(F,L,I,V)]{1,2}[E|D][(K,R)]{1,2}', '[(F,L,I,V)]{1,2}[E|D][(K,R)]{1,2}Q',\
#          'Q[E|D][(F,L,I,V)]{1,2}[(K,R)]{1,2}', '[E|D][(F,L,I,V)]{1,2}[(K,R)]{1,2}Q',\
#          'Q[E|D][(K,R)]{1,2}[(F,L,I,V)]{1,2}', '[E|D][(K,R)]{1,2}[(F,L,I,V)]{1,2}Q']
# CMA_Substrate_Dict = read_old_fasta('KillifishResources/longest-Protein-seqs_Nfu-ncbi100_andMito.fasta', motif = motif)
# XP_list = [i[-14:] for i in CMA_Substrate_Dict.keys()]
# np.savetxt('GeneSets/Killifish/CMA_Substrate.csv', XP_list, delimiter =',', fmt='%s')

# ########## Section #2: Given list of XP, output raw MS values ##########
# # df_symbol = pd.read_csv('KillifishResources/nfur_updated_symbols_20161013.csv')
# # Find all related proteins in the fish proteome
# geneset_cateogry = 'CMA_Substrate'
# geneset_fnames = ['CMA_Substrate_Killifish.csv']
# # sample_type can be WCL, or AGG, or Prop
# sample_type = 'Prop' #AGG, WCL, Prop
# pairs = [('Young','Old'),('Young','Tert'),('Old','Tert')]
# fpath = 'GeneSets/Killifish/Autophagy'
# foutpath = 'GeneSets/MS/%s_%s'%(geneset_cateogry, sample_type)
# if not os.path.exists(foutpath): os.makedirs(foutpath)
# all_geneset = []
# geneset_list = []
# df_input = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
# for geneset_fname in geneset_fnames:
#     with open(os.path.join(fpath, geneset_fname), 'r') as f: geneset = f.read()
#     geneset = geneset.replace('\r','\n')
#     geneset = geneset.replace('\n\n','\n')
#     geneset = list(set(geneset.split('\n')))
#     all_geneset += geneset
#     geneset_list.append(geneset)
# for pair in pairs:
#     age_pair = '%sv%s'%(pair[1][0], pair[0][0])
#     for geneset, set_name, i in zip(geneset_list, geneset_fnames, range(len(geneset_list))):
#         #### generate table that's convenient for heatmap or clustering ####
#         if sample_type == 'Prop':
#             df_data = df_input[(df_input['N. furzeri Protein Id'].isin(geneset)) & (df_input['Sample.Type']=='AGG')]
# #            fc = 'Delta_Propensity_%s-%s'%(pair[1], pair[0])
# #            pval = 'Propensity_%s-%s_pval'%(pair[1], pair[0])
#             fc = '%s_prop_logFC'%age_pair
#             pval = '%s_prop_pval'%age_pair
#         else:
#             df_data = df_input[(df_input['N. furzeri Protein Id'].isin(geneset)) & (df_input['Sample.Type']==sample_type)]
#             fc = '%s_logFC'%age_pair
#             pval = '%s_pval'%age_pair
#         if (df_data[fc].dropna(how='all').empty) == False:
#             df_data = new_rank(df_data, fc = fc, pval = pval)
#             df_out = df_data.pivot_table(index = 'Human', columns = 'Tissue', values = [fc, pval, 'rankstats'])
#             df_out.insert(0, 'classification', set_name.replace('_Killifish.csv','').replace('_All',''))
#             df_out_flat = df_data[['Human', 'Tissue', fc, pval, 'rankstats']]
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
# ############# End of Section 2 ##############


# ########## Section #3: Given list of XP, output raw MS values ##########
# df_symbol = pd.read_csv('KillifishResources/nfur_updated_symbols_20161013.csv')
# Find all related proteins in the fish proteome
fpath = 'GeneSets/Killifish/Autophagy'
geneset_cateogry = 'CMA_Substrate'
geneset_fname = 'CMA_Substrate_Killifish.csv'
geneset = open(os.path.join(fpath, geneset_fname), 'r').read().split('\n')
# sample_type can be WCL, or AGG, or Prop
sample_type = 'AGG' #AGG, WCL, Prop
pairs = [('Young','Old'),('Young','Tert'),('Old','Tert')]
tissue_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']#,'Testis']
foutpath = 'GeneSets/MS/%s_%s'%(geneset_cateogry, sample_type)
if not os.path.exists(foutpath): os.makedirs(foutpath)
df_input = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
df_input['CMA_Sub']=df_input['N. furzeri Protein Id'].isin(geneset)
ttest_out = np.zeros((len(pairs), len(tissue_order)))
for p_i, pair in enumerate(pairs):
    age_pair = '%sv%s'%(pair[1][0], pair[0][0])
    if sample_type == 'Prop':
        fc = '%s_prop_logFC'%age_pair
        df_plot = df_input[df_input['Sample.Type']=='AGG']
    else:
        fc = '%s_logFC'%age_pair
        df_plot = df_input[df_input['Sample.Type']==sample_type]
    for t_i, tissue in enumerate(tissue_order):
        df_tissue = df_plot[(df_plot['Tissue']==tissue)&(df_plot[fc].notna())]
        tissue_ttest = ttest_ind(*df_tissue.groupby('CMA_Sub')[fc].apply(lambda x: list(x)), equal_var = False)
        ttest_out[p_i,t_i] = tissue_ttest[1]
    fig, ax = plt.subplots(figsize=(15,10))
    sns.boxplot(data = df_plot, x = 'Tissue', y = fc, hue = 'CMA_Sub',\
                order = tissue_order, ax = ax, palette="Set3", fliersize = 0)
    sns.swarmplot(data = df_plot, x = 'Tissue', y = fc, hue = 'CMA_Sub',\
                  order = tissue_order, ax = ax, dodge=True, size = 1, palette="Set3")
    fig.savefig(os.path.join(foutpath, 'CMA_Substrate_%s_%s.pdf'%(sample_type, age_pair)))
df_ttest = pd.DataFrame(data = ttest_out, index = pairs, columns = tissue_order)
df_ttest.to_csv(os.path.join(foutpath, 'CMA_Substrate_%s_ttestp.csv'%(sample_type)))

# # ########## Section #4: Compare brain TL and AGG substrate differences ##########
# # df_symbol = pd.read_csv('KillifishResources/nfur_updated_symbols_20161013.csv')
# # Find all related proteins in the fish proteome
# fpath = 'GeneSets/Killifish/Autophagy'
# geneset_cateogry = 'CMA_Substrate'
# geneset_fname = 'CMA_Substrate_Killifish.csv'
# geneset = open(os.path.join(fpath, geneset_fname), 'r').read().split('\n')
# pairs = [('Young','Old'),('Young','Tert'),('Old','Tert')]
# tissue = 'Brain' #_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin','Testis']
# foutpath = 'GeneSets/MS/%s_TL'%(geneset_cateogry)
# if not os.path.exists(foutpath): os.makedirs(foutpath)
# df_input = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
# df_input['CMA_Sub']=df_input['N. furzeri Protein Id'].isin(geneset)
# for p_i, pair in enumerate(pairs):
#     age_pair = '%sv%s'%(pair[1][0], pair[0][0])
#     df_plot = df_input[df_input['Tissue']==tissue]
#     fc = '%s_logFC'%age_pair
#     fig, ax = plt.subplots(figsize=(6,8))
#     sns.boxplot(data = df_plot, x = 'Sample.Type', y = fc, hue = 'CMA_Sub', order = ['TL','AGG'], ax = ax,\
#                 fliersize = 0, palette="Set3")
#     sns.swarmplot(data = df_plot, x = 'Sample.Type', y = fc, hue = 'CMA_Sub', order = ['TL','AGG'], ax = ax, \
#                   dodge=True, size = 1.5, palette='Set3')
#     fig.savefig(os.path.join(foutpath, 'CMA_Substrate_%s_%s.pdf'%(tissue, age_pair)))