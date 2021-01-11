#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 07:54:08 2019

this script is used to clean up the enrichment analysis table

@author: yiwenchen
"""

import os
import pandas as pd
import numpy as np

# this filter_core function will eliminate terms where the fil_coln has identical entries or is a subset of another fil_coln term
def filter_core(df, sort_coln='NES', fil_coln = 'core_enrichment', drop_dup = True, drop_subset = True):
    if sort_coln == 'NES':
        df['NES_abs'] = df['NES'].abs()
        sort_coln = 'NES_abs'
    # fil_coln contains the list of genes    
    df_sort = df.sort_values([fil_coln, sort_coln], axis = 0, ascending = False)
    # get rid of identical entries
    if drop_dup == True:
        df_filter = df_sort.drop_duplicates(subset=[fil_coln], keep = 'first')
        df_filter = df_filter.dropna(subset=[fil_coln], axis = 0)
        df_filter.reset_index(inplace=True, drop=True)
    else:
        df_filter = df_sort
    # get rid of subset of entries
    if drop_subset == True:
        drop_row = []
        df_filter.loc[:,'GeneCounts'] = df_filter[fil_coln].apply(lambda x: len(x.split('/')))  #############################
        df_filter = df_filter.sort_values(['GeneCounts'], axis = 0, ascending=True)
        df_filter.reset_index(inplace=True, drop=True)
        for index, row in df_filter.iterrows():
            query_list = row[fil_coln].split('/')
            if index < (len(df_filter.index) - 1):
                for index_rest, row_rest in df_filter.loc[index+1:].iterrows():
                    current_list = row_rest[fil_coln].split('/')
                    if set(query_list).issubset(set(current_list)):
                        drop_row.append(index)
                        break
        df_filter.drop(df_filter.index[drop_row],inplace=True)
    return df_filter


def table_cleanup(fpath, fname, ftype, foutpath, onto_fil = ['CC','MF'], sort_coln = 'NES', \
                  fil_coln='core_enrichment', pvalcutoff = 1, drop_dup = True, drop_subset = True):
    df_input = pd.read_csv(os.path.join(fpath, fname))
    df_input = df_input[df_input['pvalue']<=pvalcutoff]
    if ftype == 'GSEA':
        df_input = df_input[df_input['ONTOLOGY'].isin(onto_fil)]
    df_sort = filter_core(df_input, sort_coln = sort_coln, fil_coln = fil_coln, drop_dup = drop_dup, drop_subset = drop_subset)
    return df_sort
 

## note this pval_filter will keep any term where at least one of the tissue has a significant enrichment
def pval_filter(df, tissues, pvalcutoff = 0.05):
    df.loc[:,'min_pval'] = df[['pvalue_%s'%tissue for tissue in tissues]].min(axis=1, skipna=True)
    df_new = df[df['min_pval']<=pvalcutoff]
    df_new.drop(columns=['min_pval'], inplace=True)
    return df_new
 
      
# function used to apply to each tissue and count number of terms that has either significant positive or negative enrichment
def count_pval_dir(row, tissues, pvalcutoff = 0.05):
    pos_counts = 0
    neg_counts = 0
    for tissue in tissues:
        if (row['pvalue_%s'%tissue] <= pvalcutoff):
            if (row['NES_%s'%tissue] > 0):
                pos_counts += 1
            else:
                neg_counts +=1
    return pd.Series([pos_counts, neg_counts], index = ['possigpval_count','negsigpval_count'])


def count_pvalcut(df, tissues, pvalcutoff = 0.05):
    df.loc[:,'sigpval_count'] = df[df[['pvalue_%s'%tissue for tissue in tissues]]<=pvalcutoff].count(axis=1)
    df[['possigpval_count','negsigpval_count']] = df.apply(lambda row: count_pval_dir(row, tissues, pvalcutoff = pvalcutoff), axis = 1)
    return df


def inclusion(df, nbest, nsig, tissues = '', direction = 'both', cm_pfilter=0.05):
    if tissues =='': tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    if nbest > 0:
        for i, tissue in enumerate(tissues):
            df_pos = df[df['pvalue_%s'%tissue]<=cm_pfilter].nlargest(nbest, columns = ['NES_%s'%tissue])
            df_pos = df_pos[(df_pos['NES_%s'%tissue]>=0) & (df_pos['pvalue_%s'%tissue]<=cm_pfilter)]
            df_neg = df[df['pvalue_%s'%tissue]<=cm_pfilter].nsmallest(nbest, columns = ['NES_%s'%tissue])
            df_neg = df_neg[(df_neg['NES_%s'%tissue]<=0) & (df_neg['pvalue_%s'%tissue]<=cm_pfilter)]
            if i == 0:
                if direction.lower() == 'both':
                    df_inclusion = df_pos.append(df_neg,ignore_index=True)
                elif direction.lower() == 'up':
                    df_inclusion = df_pos.copy(deep =True)
                elif direction.lower() == 'down':
                    df_inclusion = df_neg.copy(deep =True)
            else:
                if (direction.lower() == 'up') or (direction.lower() == 'both'):
                    df_inclusion = df_inclusion.append(df_pos,ignore_index=True)
                if (direction.lower() == 'down') or (direction.lower() == 'both'):
                    df_inclusion = df_inclusion.append(df_neg,ignore_index=True)
    else:
        df_inclusion = pd.DataFrame()
    # nsig = minimum number of term that shared significant enrichment to be visualized
    if nsig > 0:
        if direction == 'both':
            df_inclusion = df_inclusion.append(df[df['sigpval_count']>=nsig], ignore_index=True)
        elif direction.lower() == 'up':
            df_inclusion = df_inclusion.append(df[df['possigpval_count']>=nsig], ignore_index=True)
        elif direction.lower() == 'down':
            df_inclusion = df_inclusion.append(df[df['negsigpval_count']>=nsig], ignore_index=True)
    if df_inclusion.empty == False:
        df_inclusion.drop_duplicates(keep = 'last', inplace=True)
    return df_inclusion


def f_suffix(metric, ftype, age_comp, subtype = ''):
    if ftype == 'GSEA':
        fsuffix = '%s_GO%s_%s.csv'%(metric, ftype, age_comp)
    elif subtype == '':
        fsuffix = '%s_%s_%s.csv'%(metric, ftype, age_comp)
    else:
        fsuffix = '%s_%s-%s_%s.csv'%(metric, ftype, subtype, age_comp)
    return fsuffix 


# compute a rank stats that incorporate pval as well as the enrichment score, 
# after transforming the two metric of interest should change monotonically in one direction
# a small rank stats (negative) means that the geneset is enriched among the list that show decreased expression
# a large positive rank stats means that the geneset is enriched among the list that show increased expression    
def rank_stats(df, pval, fc, stat_coln='rank_stats'):
#    df[stat_coln] = -np.log10(df[pval])*df[fc]
    df.loc[:,stat_coln] = -np.log10(df[pval])*df[fc]
    return df


# this function is used to combine all the separate GSEA results for particular database (i.e. GO, MSigDb, KEGG) 
# from mutliple tissues into one
def combine_table(fpath, fsuffix, foutpath, arrange_by, ftype, tissues='', sort_coln = 'NES', fil_coln='core_enrichment',\
                  cm_pfilter = 0.05, onto_fil = ['CC','MF'], drop_dup = False, drop_subset = False, outfname = None):
    if tissues =='': tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    if ('t' in fsuffix.lower()) & ('Testis' in tissues): tissues.remove('Testis')
    id_coln = ['ID',	 'Description']
    if ftype == 'GSEA': id_coln.append('ONTOLOGY')
    combine_coln = id_coln + ['NES', 'pvalue']
    if arrange_by not in combine_coln: combine_coln.append(arrange_by)
    for i, tissue in enumerate(tissues):
        fname = '%s_%s'%(tissue,fsuffix)
        df_tissue = pd.read_csv(os.path.join(fpath, fname))
        if ftype == 'GSEA':
            df_tissue = df_tissue[df_tissue['ONTOLOGY'].isin(onto_fil)]
        # this filter_core generally will eliminate terms where the fil_coln has identical entries or is a subset of another fil_coln term
        # but here i want to preserve all terms and filter later based on the description (i.e. similar geneset in one 
        # tissue for a TERM that would be otherwise eliminated might be the only term that showed up for another tissue)
        df_tissue = filter_core(df_tissue, sort_coln='NES', fil_coln = 'core_enrichment', drop_dup = drop_dup, drop_subset = drop_subset)
        df_tissue_add = df_tissue[combine_coln]
        # add the NES and pvalue term for each tissue
        df_tissue_add.columns = ['%s_%s'%(col, tissue) if col in ['NES', 'pvalue'] else col for col in df_tissue_add.columns]
        # the rank stats which is simply defined as -np.log10(pval)*(NES)
        df_tissue_add = rank_stats(df_tissue_add, 'pvalue_%s'%tissue, 'NES_%s'%tissue, stat_coln='rank_stats_%s'%tissue)
        if i == 0:
            df_combine_nofil = df_tissue_add  
        else:
            df_combine_nofil = pd.merge(df_combine_nofil, df_tissue_add, how = 'outer', on= id_coln)
    if outfname == None:
        df_combine_nofil.to_csv(os.path.join(foutpath, 'All_Raw_%s'%fsuffix), index = False)
    else:
        df_combine_nofil.to_csv(os.path.join(foutpath, outfname), index = False)
    return df_combine_nofil


def dtype_summary(df_raw, age_comp, direction, df_nodupsig = None, cm_pfilter=0.05, nbest = None, nsig = None, tissues = None):  
    if tissues == None: tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    if 'T' in age_comp: tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']  
    if nbest == None: nbest = 3 # number of term that are either significantly enrichment in top (increase with age) or bottom (decrease with age)
    if nsig == None: nsig = 3 # minimum number of term that shared significant enrichment to be visualized
    # fitler by number of best terms for each tissues, and minimum number of shared significant terms across tissues, add direction information here if needed
    # keep only terms where there is at least one significant pvalue
    df_raw_fil = pval_filter(df_raw, tissues, pvalcutoff = cm_pfilter)
    # figure out how many significant pvalues is the specifc GO term associated with (separate positve vs. negative enrichment score) 
    df_raw_fil = count_pvalcut(df_raw_fil, tissues, pvalcutoff = cm_pfilter)
    if df_nodupsig is not None:
        df_nodupsig_fil = pval_filter(df_nodupsig, tissues, pvalcutoff = cm_pfilter)
        df_nodupsig_fil = count_pvalcut(df_nodupsig_fil, tissues, pvalcutoff = cm_pfilter)
    # subslice the dataframe to include nbest terms and nsig shared among tissues, specify direction of changes when needed
    # for function inclusion the input variables are df, nbest, nsig in order
    if nsig != 0:
        df_raw_summary = inclusion(df_raw_fil, 0, nsig, tissues = tissues, direction = direction, cm_pfilter = cm_pfilter)
    if nbest != 0:
        # for the shared terms, to capture as much overlap, instead we use the raw table and disregard individual tissue duplicates and subset terms
        if df_nodupsig is None:
            df_nodupsig_summary = inclusion(df_raw_fil, nbest, 0, tissues = tissues, direction = direction, cm_pfilter = cm_pfilter) 
        else:
        # alternatively, one can feed a different dataframe where duplicates and subset terms are dropped
            df_nodupsig_summary = inclusion(df_nodupsig_fil, nbest, 0, tissues = tissues, direction = direction, cm_pfilter = cm_pfilter)    
    if (nbest > 0) and (nsig > 0):
        df_summary = pd.concat([df_nodupsig_summary, df_raw_summary], sort = False, ignore_index = True)
        return df_summary
    elif (nbest > 0) and (nsig == 0):
        return df_nodupsig_summary
    elif (nbest == 0) and (nsig > 0):
        return df_raw_summary


###### Section 1: Try combine all the separate GSEA results for a particular database first ###### 
#### then sort and filter based on number of significant terms ####
##### part 1.
# age_comp = 'OvY' # 'TvO', 'TvY'
# cm_pfilter = 0.05
# nbest = 3 # number of term that are either significantly enrichment in top (increase with age) or bottom (decrease with age)
# nsig = 0 # minimum number of term that shared significant enrichment to be visualized
# #directions = ['both']
# #directions = ['up','down', 'both']
# for direction in ['down','up']:
#     foutpath = "GSEA/Cleanedup/%s"%age_comp
#     if not os.path.exists(foutpath): os.makedirs(foutpath)
#     arrange_by = 'NES' 
#     #### go through different table styles ####
#     fstype_dict = {'GSEA':[''], 'KEGG':['', 'Modules'], 'MSigDb': ['GO-ALL', 'GOCC', 'HALLMARKS'], 'DO':['']}
#     ftypes = ['GSEA', 'KEGG', 'MSigDb', 'DO']
#     #### section 2.1: use raw data (drop_dup = False, drop_subset = False) ####
#     for ftype in ftypes:
#         for metric  in ['AGG', 'TL', 'Prop']:
#             subtypes = fstype_dict[ftype]
#             for subtype in subtypes:
#                 fpath = "GSEA/Results/%s/%s/%s"%(metric, age_comp, ftype)
#                 fsuffix = f_suffix(metric, ftype, age_comp, subtype = subtype)
#                 raw_fname = os.path.join(foutpath, 'All_Raw_%s'%fsuffix)
#                 nodupsig_fname = 'All_nodupsub_%s'%fsuffix
#                 # read the general combined table first, if it is not available then create that
#                 # note that raw_fname = 'All_Raw_%s'%fsuffix are files without dropping duplicates and terms with genes that could be a subset of another
#                 if os.path.isfile(raw_fname): 
#                     df_raw = pd.read_csv(raw_fname)
#                 else:
#                     df_raw = combine_table(fpath, fsuffix, foutpath, arrange_by, ftype, onto_fil = ['MF'], drop_dup = False, drop_subset = False)
#                 # note that 'All_nodupsub_%s'%fsuffix are files where duplicates and terms with genes that could be a subset of another dropped  
#                 df_vis = dtype_summary(df_raw, age_comp, direction, cm_pfilter = cm_pfilter, nbest = nbest, nsig = nsig)
#                 df_vis.to_csv(os.path.join(foutpath, 'All_%s_nbest%s_nsig%s_%s_fromraw.csv'%(fsuffix.replace('.csv',''), nbest, nsig, direction.lower())),index = False)
# #### part 2. use the following code if you want to drop gene set that is part of another gene set ####
# #### this section will take a while to run especially foor MSigDb results ####
# age_comp = 'OvY' # 'TvO', 'TvY'
# cm_pfilter = 0.05
# nbest = 3 # number of term that are either significantly enrichment in top (increase with age) or bottom (decrease with age)
# nsig = 0 # minimum number of term that shared significant enrichment to be visualized
# #directions = ['both']
# #directions = ['up','down', 'both']
# for direction in ['up']:
#     foutpath = "GSEA/Cleanedup/%s"%age_comp
#     arrange_by = 'NES' 
#     #### go through different table styles ####
#     fstype_dict = {'GSEA':[''], 'KEGG':['', 'Modules'], 'MSigDb': ['GO-ALL', 'GOCC', 'HALLMARKS'], 'DO':['']}
#     ftypes = ['GSEA', 'KEGG', 'MSigDb', 'DO']
#     #### section 2.1: use raw data (drop_dup = False, drop_subset = False) ####
#     for ftype in ftypes:
#         for metric  in ['AGG', 'TL', 'Prop']:
#             subtypes = fstype_dict[ftype]
#             for subtype in subtypes:
#                 fpath = "GSEA/Results/%s/%s/%s"%(metric, age_comp, ftype)
#                 fsuffix = f_suffix(metric, ftype, age_comp, subtype = subtype)
#                 raw_fname = os.path.join(foutpath, 'All_Raw_%s'%fsuffix)
#                 nodupsig_fname = 'All_nodupsub_%s'%fsuffix
#                 # read the general combined table first, if it is not available then create that
#                 # note that raw_fname = 'All_Raw_%s'%fsuffix are files without dropping duplicates and terms with genes that could be a subset of another
#                 if os.path.isfile(raw_fname): 
#                     df_raw = pd.read_csv(raw_fname)
#                 else:
#                     df_raw = combine_table(fpath, fsuffix, foutpath, arrange_by, ftype, onto_fil = ['MF'], drop_dup = False, drop_subset = False)
#                 # note that 'All_nodupsub_%s'%fsuffix are files where duplicates and terms with genes that could be a subset of another dropped  
#                 if os.path.isfile(os.path.join(foutpath, nodupsig_fname)):
#                     df_nodupsig = pd.read_csv(os.path.join(foutpath, nodupsig_fname))
#                 else:
#                     df_nodupsig = combine_table(fpath, fsuffix, foutpath, arrange_by, ftype, onto_fil = ['MF'], drop_dup = True, drop_subset = True, outfname = nodupsig_fname)
#                 df_vis = dtype_summary(df_raw, age_comp, direction, df_nodupsig = df_nodupsig, cm_pfilter = cm_pfilter, nbest = nbest, nsig = nsig)
#                 df_vis.to_csv(os.path.join(foutpath, 'All_%s_nbest%s_nsig%s_%s.csv'%(fsuffix.replace('.csv',''), nbest, nsig, direction.lower())),index = False)              
### End of Section 1 ######