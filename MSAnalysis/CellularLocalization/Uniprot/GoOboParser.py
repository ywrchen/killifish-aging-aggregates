#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:00:09 2019

# this is the script to 
1. parse GO ontology .obo file and make it into a table go_obo_table.csv
2. read the goa_human.gaf file
3. merge the goa_human.gaf file with the .obo file to get the GO name explicitly
4. add on the ensembl protein id to match with the corresponding uniprot id in goa_human, 
    the resulting file is 'goa_human_ensembl.csv'
5. filter out the specific go term for all killifish proteins based on its human homolog
    the resulting file is 'goa_killifish.csv'
    
# for details about the .obo file, see http://geneontology.org/docs/download-ontology/

@author: yiwenchen
"""

import os
import pandas as pd

def read_obo(fpath, fname):
    f = open(os.path.join(fpath, fname), 'r')
    entries = f.read().split('\n\n')
    data_dict = dict()
    all_keys = []
    header = entries[0]
    for i, go_i in enumerate(entries[1:]):
        go_row = go_i.split('\n')[1:]
        row_keys = []
        if '' in go_row: go_row.remove('')
        for j in go_row:
            j_colon_index = j.index(': ')
            j_key = j[:j_colon_index]
            j_entry = j[j_colon_index+2:]
            if j_key not in all_keys:
                all_keys.append(j_key)
                row_keys.append(j_key)
                data_dict[j_key] = ['']*i + [j_entry]
            elif j_key not in row_keys:
                row_keys.append(j_key) 
                data_dict[j_key].append(j_entry)
            elif j_key in row_keys:
                data_dict[j_key][i] += '|%s'%(j_entry)
        for empty_key in (set(all_keys) - set(row_keys)):
            data_dict[empty_key].append('')
    df_obo = pd.DataFrame(data_dict)
    return df_obo[all_keys], header


def filter_header(fpath, fname):
    f = open(os.path.join(fpath, fname), 'r')
    i = 0
    while '!' in f.readline():
        i += 1
    return i


def DB_clean(df_query, df_filter, DB):
    df_filter['Assigned By'] = DB
    # merge the filter to query to apply replacement entry easily
    df_clean = pd.merge(df_query, df_filter, how = 'left', left_on = ['Assigned By', 'name'],\
                        right_on = ['Assigned By', 'name'])
    df_clean.loc[df_clean['Keep']== False, 'name'] = df_clean['Replace']
    # remove used columns
    df_clean.drop(columns = ['Keep', 'Replace'], inplace = True)
    # get rid of duplicates seen with the same "Assigned By" and "name" entry
    df_clean.drop_duplicates(subset = ['Kf Id', 'Assigned By', 'name', 'namespace'], keep = 'first', inplace=True)
    df_clean.dropna(how = 'any', subset = ['name'], inplace=True)
    # clean up the query further to use only Uniprot and Reactome entries and keep only cellular_components term
    df_out = df_clean[(df_clean['namespace'] == 'cellular_component') & (df_clean['Assigned By'].isin(['UniProt', 'Reactome']))]
    return df_out


def clean_DB_term(row, df_all):
    Kf_Id = row['Kf Id']
    all_terms = df_all['name'][df_all['Kf Id'] == Kf_Id]
    term_i = row['name']
    if term_i in ['nucleus', 'cytoplasm']:
        if len(all_terms) >1:
            all_terms = all_terms.tolist()
            if ('nucleus' in all_terms) and ('nucleoplasm' in all_terms) and (term_i == 'nucleus'):
                row['name'] = 'nucleoplasm'
            if ('cytoplasm' in all_terms) and ('cytosol' in all_terms) and (term_i == 'cytoplasm'):
                row['name'] = 'cytosol'
    return row
    

def CM_clean(df_query):
    # remove duplicate entries from UniProt and Reactome if they are the same
    df_clean = df_query.drop_duplicates(subset = ['Kf Id', 'name', 'namespace'], keep = 'first')
    df_ref = df_clean.copy() 
    # remove higher order terms from lower order terms (i.e. protein mapped to nucleus and nucleoplasm, only keep nucleoplasm)
    df_out = df_clean.apply(clean_DB_term, axis = 1, args = (df_ref,))
    df_out.drop_duplicates(subset = ['Kf Id', 'name', 'namespace'], keep = 'first', inplace=True)    
    return df_out

############ SECTION 1: RUN THIS TO RE-FORMAT OBO FILE TO TABLE ########### 
fpath = ''
fname = 'go.obo'
df_obo, header = read_obo(fpath, fname)
df_obo.to_csv('go_obo_table.csv', index = False)
############################# END OF SECTION 1 ############################


######### SECTION 2: RUN THIS TO MERGE GAF, DOBO, and Ensembl Info ######## 
df_obo = pd.read_csv('go_obo_table.csv')
# read the goa_human.gaf file
fpath = ''
f_gaf = 'goa_human.gaf'
coln_gaf = ['DB', 'Uniptot ID', 'DB Object Symbol', 'Qualifier', \
            'GO ID', 'DB:Reference', 'Evidence Code', 'With (or) From',\
            'Aspect', 'DB Object Name', 'DB OBject Synonym', 'DB Object Type', \
            'Taxon', 'Date', 'Assigned By', 'Annotation Extension', 'Gene Product Form ID']
human_gaf_header_end = filter_header(fpath, f_gaf)
df_human_gaf = pd.read_csv(f_gaf, delimiter = '\t', \
                          skiprows = human_gaf_header_end, header = None, names = coln_gaf)
# merge goa_human.gaf with go_obo_table
df_human_go = df_human_gaf.merge(df_obo, how = 'left', left_on = 'GO ID', right_on = 'id')
df_human_go.drop('id', axis=1, inplace=True)
# merge the ensembl id with uniprot ide to the go_table
df_ensembl_uniprot = pd.read_csv('ensembl-uniprot_human.tab', delimiter = '\t', header = None, \
                                names = ['Uniptot ID', 'Human Ensembl Protein Id'])
df_human_go = df_human_go.merge(df_ensembl_uniprot, how = 'left', left_on = 'Uniptot ID', right_on = 'Uniptot ID')
df_human_go.to_csv('goa_human_ensembl.csv', index = False)
############################# END OF SECTION 2 #############################


######### SECTION 3: RUN THIS TO MERGE KILLIFISH HUMAN HOMOLOGS TO GO #######
df_human_go = pd.read_csv('goa_human_ensembl.csv')
## the human ortholog ensembl id were obtained from BIOMART and merged
df_kf = pd.read_csv('Killifish_human_orthologs.csv') 
df_kf_go = df_kf.merge(df_human_go, how = 'left', left_on = 'Human Ensembl Protein Id', right_on = 'Human Ensembl Protein Id')
df_kf_go.to_csv('killifish_human_go.csv', index = False)
## merge the human go to ensembl protein name and uniprot id
############################## END OF SECTION 3 #############################


########## SECTION 4: GET ALL UNIPROT CELLULAR LOCALIZATION TERMS ###########
df_kf_go = pd.read_csv('killifish_human_go.csv')
df_uniprot_count = df_kf_go[(df_kf_go['Assigned By'] == 'UniProt') & (df_kf_go['namespace'] == 'cellular_component')].groupby('name').count()
df_uniprot_count.to_csv('uniprot_counts.csv')
############################## END OF SECTION 4 #############################


####### SECTION 5: CONSOLIDATE UNIPROT AND REACTOME LOCALIZATION TERMS ######
df_kf_go = pd.read_csv('killifish_human_go.csv')
# uniprot_CM_cleanup.csv and reactome_CM_cleanup are manually curated cellular compartment terms to minimize redundancy.
df_uniprot_filter = pd.read_csv('uniprot_CM_cleanup.csv')
df_reactome_filter = pd.read_csv('reactome_CM_cleanup.csv')
df_clean = DB_clean(df_kf_go, df_uniprot_filter, 'UniProt')  
df_clean = DB_clean(df_clean, df_reactome_filter, 'Reactome')
df_kf_localization = CM_clean(df_clean) 
df_kf_localization.to_csv('killifish_human_CM.csv', index = False) 
############################## END OF SECTION 5 #############################