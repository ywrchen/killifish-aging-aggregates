#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 16:42:03 2021

@author: yiwenchen

this is to create genelist that serves as input for GSEA analysis
"""

import pandas as pd
import os
import numpy as np


# read the input file
df_all = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
metrics = ['TL', 'AGG', 'Prop'] 
age_pairs =  ['OvY', 'TvO', 'TvY']
# loop through metrics, then age_pair, then tiissues to get the rank-stats
if not os.path.exists('GSEA/GeneSets'): os.makedirs('GSEA/GeneSets')
for metric in metrics:
    if not os.path.exists('GSEA/GeneSets/%s'%metric): os.makedirs('GSEA/GeneSets/%s'%metric)
    if metric == 'TL':
        df_metric = df_all[df_all['Sample.Type']==metric]
    else:
        df_metric = df_all[df_all['Sample.Type']=='AGG']
    for age_pair in age_pairs:
        if 'T' in age_pair:
            tissues_order = ['Brain','Gut','Heart','Liver','Muscle','Skin']
        else:
            tissues_order = ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']   
        for tissue in tissues_order:
            id_coln = 'N. furzeri Protein Id'
            p_coln = '%s_pval'%age_pair if metric !='Prop' else '%s_prop_pval'%age_pair
            logfc_coln = '%s_logFC'%age_pair if metric !='Prop' else '%s_prop_logFC'%age_pair
            df_tissue = df_metric[[id_coln, p_coln, logfc_coln]][df_metric['Tissue']==tissue]
            rankstas = -np.log10(df_tissue[p_coln])*(df_tissue[logfc_coln])
            df_tissue = pd.concat([df_tissue[id_coln], rankstas], axis=1)
            df_tissue.columns = ['Protein', 'mlog10QvalxFC']
            df_tissue.to_csv(os.path.join('GSEA/GeneSets/%s/%s'%(metric, age_pair), '%s_%s_%s.csv'%(tissue, age_pair, metric)), index=False)