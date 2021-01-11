"""

This is the main script used to compile and analyze all the mass spectrometry data listed in MSResults/BGGLPM folder. 
There are 14 files organized by tissue origin (brain, gut, heart, liver, muscle, skin, and testis) and sample types 
(tissue lysate (TL) or aggregate (AGG)). 

Before running the script, make sure TMTUtl.py is in the same working directory

This file requires input from command line. So make sure to navigate to the current working directory and make sure
the input folder is accessible in subdirectories.

Various homolog information is mapped onto the output so make sure 'nfur_updated_symbols_20161013.csv' is in the 
working directory.

Below is an example command line input. You can chooose to analyze all files (see below) or a specific tissue.

# Command line arguments
# cd killifish-aging-aggregates/MSResults

# to Analyze all tissues
# python TissueTMT.py --FileName BGLPM/*_TargetProtein.txt

# to Analyze brain only
# python TissueTMT.py --FileName BGLPM/Brain_AGG_PMOnly_TargetProtein.txt BGLPM/Brain_TL_PMOnly_TargetProtein.txt

"""

import argparse
import itertools
import matplotlib
from pylab import rcParams
rcParams.update(matplotlib.rcParamsDefault)
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import os
import os.path
import pandas as pd
import time

import TMTUtl as tmt_utl


#rcParams['figure.figsize'] = 20, 20

################################# PARSE INPUT AND STORE THEM TO INITIALIZE VARIOUS PARAMETERS ##############################
# Read the command line arguments
parser = argparse.ArgumentParser(description='Process input files and analyze the data')
#parser.add_argument('FileName', nargs = '+', default = os.listdir(os.getcwd()),
                    #help = 'The name of the protein files (has to contain tissue type and sample type (AGG or TL), best if in the format of  ProteomeType-TissueType, can use * to simplify inputs,')
#parser.add_argument('--FileName', nargs = '+', default = ['BGLPM/'+i for i in os.listdir('BGLPM/')],
parser.add_argument('--FileName', nargs = '+',
                    help = 'The name of the protein files (has to contain tissue type and sample type (AGG or TL), best if in the format of  ProteomeType-TissueType, can use * to simplify inputs,')
parser.add_argument('--FishName',nargs = '+', default = [],
                    help = 'the name of the fish to be analyzed')
parser.add_argument('--DropTags',nargs = '+', default = ['TL_Testis_129_C', 'TL_Testis_130_N', 'TL_Testis_130_C', 'AGG_Testis_129_C', 'AGG_Testis_130_N', 'AGG_Testis_130_C'],
                    help = 'the sample tags that should be excluded due to low signal etc. The format should be sample_type_tissue_tag')
parser.add_argument('--YoungTags', nargs = '+', default = ['128_N', '128_C', '129_N'],
                   help = 'the default mass tag for old fish is [128_N, 128_C, 129_N], change if it is different')
parser.add_argument('--OldTags', nargs = '+', default = ['126', '127_N', '127_C'],
                    help = 'the default mass tag for old fish is [126, 127_N, 127_C], change if it is different')
parser.add_argument('--TertTags', nargs = '+', default = ['129_C', '130_N', '130_C'],
                   help = 'the default mass tag for old fish is [129_C, 130_N, 130_C], change if it is different')
parser.add_argument('-PoolTag', default = '131',
                   help = 'the default mass tag for the pool channel is 131, change it is different')
parser.add_argument('--MinRatioSum', default = 600,
                   help = 'please specify the minium sum of ratio for a given channel, the default is 800')
parser.add_argument('--PSMsNormSum', default = 100000,
                   help = 'please specify the total sum of ratio for a given channel, the default is 100000')
parser.add_argument('--OutDir', default = '',
                   help = 'please specify the out directory name')
parser.add_argument('--ExcludedAgeCombs', nargs = '+', default = ['Testis-1', 'Testis-2'],
                   help = 'please specify the age combination that does not exit or shouldnt be analyzed')                   

########################### User specified input ###################################################### 
# normalization methods: 1) by sum; 2) by smallest PSMs; 3) by median PSMs; 4)by smallest PSMs/ProtLen. All on a sample/channel basis
norm_method = 1
exception = ['Accession', 'Description', '# AAs', 'MW [kDa]', 'calc. pI']
two_age_combination = [('Young','Old'), ('Young','TERT'), ('Old', 'TERT')]

# Parse the command line arguments into variables
file_names = parser.parse_args().FileName
drop_tags = parser.parse_args().DropTags
excluded_tags = [': '.join(i.split('-')) for i in drop_tags]
young_tags = parser.parse_args().YoungTags
old_tags = parser.parse_args().OldTags
tert_tags = parser.parse_args().TertTags
pool_tag = parser.parse_args().PoolTag
min_signal_sum = int(parser.parse_args().MinRatioSum)
PSMsNorm_sum = int(parser.parse_args().PSMsNormSum)
figdir = parser.parse_args().OutDir
excluded_age_combs = [(i.split('-')[0], int(i.split('-')[1])) for i in parser.parse_args().ExcludedAgeCombs]


# create output directory based on input files and parameters
if figdir == '' and '/' in file_names[0]:
    figdir = '%sOut'%file_names[0].split('/')[0]
if norm_method == 1:
    figdir = '%s_NormBySum'%figdir
elif norm_method == 2:
    figdir = '%s_NormByMin'%figdir
elif norm_method == 3:
    figdir = '%s_NormByMedian'%figdir
if not os.path.exists(figdir):
    os.makedirs(figdir)

# Create additional variable when needed
tissue_types = []
sample_types = []
sample_tags = [young_tags, old_tags, tert_tags]
age_groups = ['Young', 'Old', 'TERT']
sample_tags_dict = {'Young': young_tags, 'Old': old_tags, 'TERT': tert_tags}
tag_to_age_dict = {'128_N':'Young-1', '128_C':'Young-2', '129_N':'Young-3', '126':'Old-1', '127_N':'Old-2', '127_C':'Old-3', '129_C':'TERT-1', '130_N':'TERT-2', '130_C':'TERT-3'}
sample_tags_nogroup = list(itertools.chain.from_iterable(sample_tags))


############## LOAD ALL FILES AND COMPILE ALL THE RAWDATA INTO ONE SINGLE DATAFRAME {df_all} ############
start_time = time.time()
#print min_signal_sum
# Load TL and Agg data into a dictionary of each exp dataframe, use {exception} as index, add file description (i.e. tissue type and sample type) to indicate the origin of the data
df_exps_dict, all_sample_types, all_tissue_types, drop_tags = tmt_utl.file_to_df(file_names, exception, sample_tags_nogroup, pool_tag, min_signal_sum, excluded_tags = drop_tags)
if tissue_types == []:
    tissue_types = all_tissue_types
if sample_types == []:
    sample_types = all_sample_types
df_tissues = tmt_utl.merge_agg_tl(df_exps_dict, tissue_types, exception)


file_load_time = time.time()
print ('it takes %s s to load all files into dataframe'%(file_load_time - start_time))


######################################### Process Data by Tissue ########################################
for i, tissue in enumerate(df_tissues.keys()):
    time_begin = time.time()
    ######################## DATA NORMALIATION AND PRE-PROCESSING FOR DOWNSTREAM ANALYSIS #########################
    # Calculate the total peptide for each sample and add a new column for these values
    df_tissue = df_tissues[tissue]
    df_tissue = tmt_utl.df_preprocess(df_tissue, tissue, norm_method, sample_tags_nogroup, pool_tag, tag_to_age_dict, PSMsNorm_sum, excluded_tags = drop_tags)
    time_1 = time.time()
    print ('%s: time it takes to normalize and pre-process the raw data is %s s'%(tissue, (time_1 - time_begin)))
    ########################################### DATA VISUALIZATION ################################################
    df_tissue = tmt_utl.tl_agg_enrichment(df_tissue, tissue, sample_tags_nogroup, pool_tag, tag_to_age_dict, sample_tags, age_groups, sample_types, excluded_tags = drop_tags, dir = figdir)
    time_2 = time.time()
    print ('%s: time it takes to compare agg enrichment from tl is %s s'%(tissue, (time_2 - time_1)))
    df_tissue = tmt_utl.age_comp(df_tissue, tissue, two_age_combination, tag_to_age_dict, sample_tags_dict, excluded_combs = excluded_age_combs, dir = figdir)
    df_tissue_out = tmt_utl.split_df_tissue(df_tissue, tissue, exception)
    df_tissue_out.to_csv(os.path.join(figdir, 'TMT_All_%s.csv'%tissue), index = False)
    if i == 0:
        df_all = df_tissue_out
    else:
        df_all = pd.concat([df_all, df_tissue_out], ignore_index=True, sort=False)
    time_3 = time.time()
    print ('%s: time it takes to look for aggregation propensity change is %s s'%(tissue, (time_3 - time_2)))


#### Export all the data across tissues to csv file
all_colns = ['Tissue', 'Sample.Type','N. furzeri Protein Id']
all_colns += ['%s-%s_PSMsNorm'%(age,i) for age in age_groups for i in range(1,4,1)]
all_colns += ['log2_%s-%s_PSMsNorm'%(age, i) for age in age_groups for i in range(1,4,1)]
all_colns += ['%s_%s'%(age_comp, i) for age_comp in ['OvY','TvY','TvO'] for i in ['logFC','FC','pval', 'logFC_zscore']]
all_colns += ['log2_%s-%s_prop'%(age, i) for age in age_groups for i in range(1,4,1)]
all_colns += ['%s_prop_%s'%(age_comp, i) for age_comp in ['OvY','TvY','TvO'] for i in ['logFC','FC','pval','logFC_zscore']]
all_colns += ['%s_%s_TL'%(i, age) for i in ['count', 'log2'] for age in ['Young','Old','TERT']]
all_colns += ['log2_%s%s'%(age,i) for age in ['Young','Old','TERT'] for i in ['',  '_zscore', '_zscore_TL','_prop', '_prop_zscore']]
df_all_out = df_all[all_colns]
df_all_out.columns=[i.replace('_PSMsNorm','') for i in df_all_out.columns]
# fill in protein symbol (Param's homolog table)
df_symbol = pd.read_csv('nfur_updated_symbols_20161013.csv')
df_all_out = df_symbol.iloc[:,1:].merge(df_all_out, on = 'N. furzeri Protein Id', how='right')
df_all_out.to_csv(os.path.join(figdir, 'kf_combined_prop.csv'), index = False)

finish_time = time.time()
print ('time it takes to run the entire script is %s'%(finish_time - start_time))
