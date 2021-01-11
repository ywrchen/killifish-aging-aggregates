import sys
import numpy as np
import itertools
#from itertools import chain
#from collections import Iterable
import os.path
#import copy
#import math
#import time
#import openpyxl
#import pylab
#import matplotlib
#from matplotlib import pyplot as plt
#from matplotlib.patches import Circle, Ellipse
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d import proj3d
#from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3, venn3_circles, venn3_unweighted
import scipy.stats
#import seaborn as sns
#from pylab import rcParams
#rcParams.update(matplotlib.rcParamsDefault)
#rcParams['figure.figsize'] = 20, 20
import pandas as pd


############################# VARIOUS COLUMN HEADER/VARIABLE FORMAT ####################
############### use ':' when it is for individual sample else  use '_' #################

# TL_Gut_Ratios: (126) / (131)
def Ratio_tag(tag, pool_tag):
    ratio_tag_output = 'Ratios: (%s) / (%s)' %(tag, pool_tag)
    return ratio_tag_output

# add '>ZGC:' to a number for the FASTA file line
def zgc_number(x):
    if str(x).replace('>','').isdigit():
        return '>ZGC:' + str(x).replace('>','')
    else:
        return str(x)

# get rid of common contaminant:
def drop_contaminant(x):
    if 'contaminant' in str(x):
        return np.nan
    else:
        return x

def round_value(value):
    if value == '' or value == np.nan or str(value) == 'nan':
        value = 0
    else:
        value = int(round(value))
    return value

def df_zscore(df, coln, log = 2):
    if log == 2:
        df_coln = np.log2(df[coln])
    elif log == 10:
        df_coln = np.log10(df[coln])
    elif log == 'ln':
        df_coln = np.log(df[coln])
    else:
        df_coln = df[coln]
    coln_avg = df_coln.mean()
    coln_std = df_coln.std()
    df_output = df_coln.apply(lambda x: (x - coln_avg)/coln_std)
    return df_output

def zscore(x, x_mean, x_std):
    z = (x-x_mean)/x_std
    return z

# create a function to fill each row of a datafram with 0 only when there are noNaN value in that row from other column, so
#['a','b', 'c'
#  1 , -2, na
#  2, na, na,
# na, na, na,
#]
#will be
#['a','b', 'c'
#  1 , -2, 0
#  2, 0, 0,
# na, na, na,
#] after applying this method, filldict = {'a':#replacena, 'b': #toreplacenan}
def fillna_rowsum(df, filldict = ''):
    #    df = pd.to_numeric(df, errors = 'coerce')
    if filldict == '':
        df.fillna(0)
    else:
        for index, row in df.iterrows():
            if str(df.loc[index, :].sum()) != 'nan':
                df.loc[index,:] = df.loc[index, :].fillna(filldict, inplace = False)
    return df

def all_2vs1_division(df1, tags1, df2, tags2):
    all_combination = np.array([i for i in itertools.product(tags1, tags2)])
    numerator_tags = list(all_combination[:,1])
    denominator_tags = list(all_combination[:,0])
    df_slice = pd.DataFrame(df2[numerator_tags].values/df1[denominator_tags].values, columns = range(len(all_combination)), index = df1.index)
    return df_slice


# merge AGG and TL for a tissue together and return df_tissues
def merge_agg_tl(df_exps_dict, tissue_types, exception):
    exception += ['Tissue']
    df_tissues = dict()
    for tissue_type in tissue_types:
        if ('%s-AGG'%tissue_type not in df_exps_dict.keys()) or ('%s-TL'%tissue_type not in df_exps_dict.keys()):
            sys.exit('make sure both the AGG and TL files for %s is properly read'%tissue_type)
        df_AGG = df_exps_dict['%s-AGG'%tissue_type]
        df_TL = df_exps_dict['%s-TL'%tissue_type]
        df_AGG_TL = pd.merge(left = df_TL, right = df_AGG, how = 'outer', \
                              left_on = exception, right_on = exception, suffixes = ['_TL','_AGG'])
        df_tissues[tissue_type] = df_AGG_TL
    return df_tissues


# merge df_tissue into one big df:
def split_df_tissue(df_tissue, tissue, exception):
    if tissue == 'Testis':
        ages = ['Young','Old']
    else:
        ages = ['Young','Old','TERT']
    df_tissue.rename(columns={'Description':'N. furzeri Protein Id'}, inplace=True)
    TL_count_colns = ['%s_%s_TL'%(i, age) for i in ['count', 'log2'] for age in ages]
    TL_count_colns += ['log2_%s_zscore_TL'%age for age in ages]
    df_TL = df_tissue.copy()
    df_TL.insert(1, 'Sample.Type', 'TL')
    df_TL.columns = [i.replace('_TL','') for i in df_TL.columns]
    df_TL = df_TL[[i for i in df_TL.columns if ('_AGG' not in i)]]
    df_TL = df_TL.merge(df_tissue[['N. furzeri Protein Id']+TL_count_colns], on = 'N. furzeri Protein Id', how='left')
    df_AGG = df_tissue.copy()
    df_AGG.insert(1, 'Sample.Type', 'AGG')
    df_AGG.columns = [i.replace('_AGG','') for i in df_AGG.columns]
    df_AGG = df_AGG[[i for i in df_AGG.columns if ('_TL' not in i)]]
    df_AGG = df_AGG.merge(df_tissue[['N. furzeri Protein Id']+TL_count_colns], on = 'N. furzeri Protein Id', how='left')
    df_final = pd.concat([df_TL, df_AGG], ignore_index=True, sort=False)
    # get rid of protein not detected in certain fraction
    tag = ['Young-%s_PSMsNorm'%i for i in range(1,4,1)]
    df_final.dropna(subset=tag, how='all', inplace=True)
    df_final = df_final[df_final['N. furzeri Protein Id'].str.startswith('XP_')]
    return df_final


# Add description that specify tissue type and sample type in column header for each file and return it as dataframe
def add_description(file, exception):
    sample_type, tissue, file = file_describe(file)
    if '.txt' in file:
        df = pd.read_table(file)
    elif '.csv' in file:
        df = pd.read_csv(file)
    df.insert(0, 'Tissue', tissue)
#    df['Sample.Type'] = sample_type
#    df.columns = [c if c in exception else '%s_%s_%s'%(sample_type, tissue, c) for c in df.columns]
    return df


# extract sample_type and tissue from a file name such as file = 'TMT10plex_MS3_1_9_TargetProtein-AGG-Liver.txt'
def file_describe(file):
    tissues = ['Brain', 'Liver', 'Gut', 'Heart', 'Muscle', 'Testis', 'Skin', 'Plasma']
    file = file.lower()
#    print file
    if 'agg' in file:
        sample_type = 'AGG'
    elif 'tl' in file:
        sample_type = 'TL'
    else:
        sys.exit('make sure the file you specified contain information on sample type (AGG or TL)')
    tissue_type_find = 0
    for tissue in tissues:
        if tissue.lower() in file:
            tissue_type = tissue
            tissue_type_find += 1
    if tissue_type_find != 1:
        tissue_type = file[file.index('-') + 5 : file.index('.')]
    if '/' in file:
        file_list = file.split('/')
        file_name = file_list[-1]
        dir = '/'.join(file_list[:-1])
        file = os.path.join(dir, file_name)
    return sample_type, tissue_type, file


# read each TMT experiment files and store the dataframe into df_exps_dict as an entry mapped by 'tissue-sample_type'
def file_to_df(file_names, exception, sample_tags_nogroup, pool_tag, min_signal_sum, excluded_tags = []):
    df_exps_dict = dict()
    sample_types = []
    tissue_types = []
    for file in file_names:
        sample_type, tissue, file = file_describe(file)
        if sample_type not in sample_types:
            sample_types.append(sample_type)
        if tissue not in tissue_types:
            tissue_types.append(tissue)
        df = add_description(file, exception)
        # drop contaminant entries
        df = df[df['Description']!='(Common contaminant protein)']
        # remove all unnecessary spaces to avoid overcounting entries when merging
        for exception_col in exception:
            df[exception_col] = df[exception_col] = df[exception_col].apply(lambda r: str(r).replace(' ',''))
        # Set empty channel signal to NaN (empty channel are the ones with a total signal ratio less than min_signal_sum)
        for tag in sample_tags_nogroup:
            column_tag = 'Ratios: (%s) / (%s)' %(tag, pool_tag)
#            print '%s %s has a sum of %d'%(file, column_tag, df[column_tag].sum())
            if df[column_tag].sum() < min_signal_sum:
#                print df[column_tag].sum()
#                print '%s does not have signal' %tag
                df[column_tag] = np.nan
                df['Ratio Standard Errors [%%]: (%s) / (%s)' %(tag, pool_tag)] = np.nan
                excluded_tag = '_'.join([sample_type, tissue, tag])
                if excluded_tag not in excluded_tags:
                    excluded_tags.append(excluded_tag)
        df_exps_dict['-'.join([tissue, sample_type])] = df
    return df_exps_dict, sample_types, tissue_types, excluded_tags


def ttest_func(row, pop1, pop2):
    out = scipy.stats.ttest_ind(row[pop1], row[pop2])
    return out


################################## MAIN UTILITY FUNCTIONS FOR MAIN.PY FILE #######################################


# Data normalization and pre-processing for downstream analysis
def df_preprocess(df, tissue, norm_method, sample_tags_nogroup, pool_tag, tag_to_age_dict, PSMsNorm_sum, excluded_tags = [], sample_types = ['TL', 'AGG']):
    # Calculate the total peptide for each sample and add a new column for these values
    for sample_type in sample_types:
        ratio_tags = ['Ratios: (%s) / (%s)_%s' %(tag, pool_tag, sample_type) for tag in sample_tags_nogroup]
        all_ratio_sum = df[ratio_tags].sum(axis = 1) + 1
        for tag in sample_tags_nogroup:
            ratio_tag = 'Ratios: (%s) / (%s)_%s' %(tag, pool_tag, sample_type)
            tag_label = tag_to_age_dict[tag]
            if '_'.join([sample_type, tissue, tag]) not in excluded_tags:
                PSMs_label = '%s_PSMs_%s'%(tag_label, sample_type)
                df[PSMs_label] = df[ratio_tag] * (df['# PSMs_%s'%sample_type]/all_ratio_sum)
                Norm_label = '%s_PSMsNorm_%s'%(tag_label, sample_type)
                # Normalize each channel (approach 1: normalize by sum of total signal per channel, approach 2: normalize by the min_PSMs of each channel, aroach 3: normalize to median)
                if norm_method == 1:
                    # Normalization Approach 1: normalize each channel based on sum of matched peptide counts per channel and insert a new column to indicate the normalized peptide counts scaling to an arbitrary number (PSMsNorm_sum).
                    # AGG_Gut: old-1_PSMsNorm
                    df[Norm_label] = df[PSMs_label] / (df[PSMs_label].sum()) * PSMsNorm_sum
                elif norm_method == 2:
                    # Normalization Approach 2 (by smallest PSMs for each sample)
                    df[Norm_label] = df[PSMs_label] / (df[PSMs_label].min())
                elif norm_method == 3:
                    # Normalization Approach 3 (by median PSMs for each sample)
                    df[Norm_label] = df[PSMs_label] / (df[PSMs_label].median())
                elif norm_method == 4:
                    # Normalization Approach 4 (by smallest PSMs/ProtLen for each sample)
                    PSMsOverAA_label = PSMs_label.replace('PSMs','PSMs/AA')
                    df[PSMsOverAA_label] = df[PSMs_label].divide(df['# AAs'], axis = 'index')
                    df[Norm_label] = df[PSMs_label] / (df[PSMsOverAA_label].min())
                df['log2_%s_PSMsNorm_%s'%(tag_label, sample_type)] = np.log2(df[Norm_label])
    df['# AAs'] = pd.to_numeric(df['# AAs'], errors = 'coerce')
    return df


# Analysis on how AGG enrichment is different from TL for each sample.
def tl_agg_enrichment(df, tissue, sample_tags_nogroup, pool_tag, tag_to_age_dict, sample_tags, age_groups, sample_types, excluded_tags = [], dir = ''):
    for tag in sample_tags_nogroup:
        if ('_'.join(['TL', tissue, tag]) not in excluded_tags) and ('_'.join(['AGG', tissue, tag]) not in excluded_tags):
            tag_label = tag_to_age_dict[tag]
            # AGG/TL method 1: ignore NaN_TL, for {nonNaN_TL and nonNaN_AGG}, compute AGG/TL = {AGG/TL} with no pseudo counts
            df['%s_prop'%tag_label] = pd.to_numeric(df['%s_PSMsNorm_AGG'%tag_label] / df['%s_PSMsNorm_TL'%tag_label], errors = 'coerce')
            df['log2_%s_prop'%tag_label] = np.log2(df['%s_prop'%tag_label])
    return df

# Look for age-related changes in AGGregation propensity (method 1: rank AGG/TL_old and AGG/TL_young respetively and compare the percentile differences, rank AGG_old and AGG_young and compare the percentile differences; method 2: welch t-test for AGG/TL_old and AGG/TL_young among biological replicates; method 3: t-test for AGG_old and AGG_young, TL_old and TL_young; method 4: chi-square or binomial test between AGG/TL_old and AGG/TL_young
def age_comp(df, tissue, two_age_combination, tag_to_age_dict, sample_tags_dict, excluded_combs = '', dir = ''):
    two_age_combination = [two_age_combination[i] for i in range(len(two_age_combination)) if (tissue, i) not in excluded_combs]
    for two_age in two_age_combination:
        PSMsNorm_age1 = ['%s_PSMsNorm'%tag_to_age_dict[tag] for tag in sample_tags_dict[two_age[0]]]
        PSMsNorm_age2 = ['%s_PSMsNorm'%tag_to_age_dict[tag] for tag in sample_tags_dict[two_age[1]]]
        # compare TL, AGG, and AGG/TL of old and young, compute the fold-change
        # just do old_avg/young_avg
        for sample_type in ['TL', 'AGG']:
            # calculate age-specific values
            age1_tags = ['%s_%s'%(tag, sample_type) for tag in PSMsNorm_age1]
            age2_tags = ['%s_%s'%(tag, sample_type) for tag in PSMsNorm_age2]
            df['count_%s_%s'%(two_age[0], sample_type)] = df[age1_tags].mean(axis=1)
            df['count_%s_%s'%(two_age[1], sample_type)] = df[age2_tags].mean(axis=1)
            df['log2_%s_%s'%(two_age[0], sample_type)] = np.log2(df[age1_tags]).mean(axis=1)
            df['log2_%s_%s'%(two_age[1], sample_type)] = np.log2(df[age2_tags]).mean(axis=1)
            df['log2_%s_zscore_%s'%(two_age[0], sample_type)] = df_zscore(df, 'log2_%s_%s'%(two_age[0], sample_type), log='linear')
            df['log2_%s_zscore_%s'%(two_age[1], sample_type)] = df_zscore(df, 'log2_%s_%s'%(two_age[1], sample_type), log='linear')
            # do age comparision
            age_comp = '%sv%s'%(two_age[1][0], two_age[0][0])
            # t-test for AGG_old and AGG_young, TL_old and TL_young respectively
            ttest_2vs1 = df.apply(lambda x: ttest_func(x, ['log2_%s'%tag for tag in age1_tags], ['log2_%s'%tag for tag in age2_tags]),\
                                  axis=1, result_type = 'expand')
            df['%s_pval_%s' %(age_comp, sample_type)] = ttest_2vs1[1]
            df['%s_FC_%s'%(age_comp, sample_type)] = pd.to_numeric(df[age2_tags].mean(axis=1)/df[age1_tags].mean(axis=1), errors = 'coerce')
            df['%s_logFC_%s'%(age_comp, sample_type)] = np.log2(df['%s_FC_%s'%(age_comp, sample_type)]) 
            df['%s_logFC_zscore_%s'%(age_comp, sample_type)] = df_zscore(df, '%s_logFC_%s'%(age_comp, sample_type), log='linear') 
        # calculate aggregation propensity differences
        age1_prop = ['%s_prop'%tag_to_age_dict[tag] for tag in sample_tags_dict[two_age[0]]]
        age2_prop = ['%s_prop'%tag_to_age_dict[tag] for tag in sample_tags_dict[two_age[1]]]
        df['log2_%s_prop'%two_age[0]] = df[age1_prop].mean(axis=1)
        df['log2_%s_prop'%two_age[1]] = df[age2_prop].mean(axis=1)
        df['log2_%s_prop_zscore'%two_age[0]] = df_zscore(df, 'log2_%s_prop'%two_age[0], log='linear')
        df['log2_%s_prop_zscore'%two_age[1]] = df_zscore(df, 'log2_%s_prop'%two_age[1], log='linear')
        fc = pd.to_numeric(df[age2_prop].mean(axis=1)/df[age1_prop].mean(axis=1), errors = 'coerce') 
        df['%s_prop_FC'%age_comp] = fc
        df['%s_prop_logFC'%age_comp] = np.log2(fc)
        df['%s_prop_logFC_zscore'%age_comp] = df_zscore(df, '%s_prop_logFC'%age_comp, log='linear')
        ttest_prop = df.apply(lambda x: ttest_func(x, ['log2_%s'%tag for tag in age1_prop], ['log2_%s'%tag for tag in age2_prop]),\
                              axis=1, result_type = 'expand')
        df['%s_prop_pval'%age_comp] = ttest_prop[1]
    return df
