#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 23:53:28 2020

# compare proteome data

@author: yiwenchen
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted, venn3, venn3_circles, venn3_unweighted
from matplotlib import rcParams

rcParams['pdf.fonttype'] = 42


# function used to plot venn diagram
def venn(set_list, label_list, figname, figpath):
    f_venn, ax_venn = plt.subplots(1,1)
    if len(label_list) == 2:
        venn2(set_list, label_list, ax=ax_venn)
    elif len(label_list) == 3:
        venn3(set_list, label_list, ax=ax_venn)
        # venn3 and venn3_circles take a 7-element list of subset sizes (Abc, aBc, ABc, abC, AbC, aBC, ABC), and draw a three-circle area-weighted venn diagram.
    ax_venn.set_title(figname)
    f_venn.suptitle(figname.replace('.pdf', ''))
    f_venn.savefig(os.path.join(figpath, figname))
  
    
# function used to split column into separate rows
# https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows    
def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:        
        res = res.reset_index(drop=True)
    return res


def David_age_agg(row):
    if (row['OvY_logFC']>=0) and (row['OvY_pval']<=0.05) and (row['TwoExp>1.5']=="YES"):
        return "YES"
    else:
        return "NO"
    

kf_data = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
kf_agg = kf_data[kf_data['Sample.Type']=='AGG']
kf_tl = kf_data[kf_data['Sample.Type']=='TL']

###### Section 1: this is used to infer overlap among age-associated aggregates
# Table S1D from Walther  2015, proteins with larger "log abs.LFQ AgWTD17 Corrected" \
# than "log abs.LFQ AgWTD06 Corrected" in Table S1D, in Walther 2015)
fname_worm = 'Walther2015_TS1D_logAbsLFQcorr_WT17greaterthanWT6.csv'
fpath_out = 'ProteomeComp/Walther2015_worm'
worm_data = pd.read_csv(os.path.join(fpath_out, fname_worm))
worm_data.dropna(subset =['Gene Name'], inplace=True)
query_coln = ['Tissue', 'Sample.Type', 'N. furzeri Protein Id', 'Human', 'Worm','log2_Young_zscore','log2_Old_zscore','OvY_logFC', 'OvY_pval']
df_merge = worm_data.merge(kf_agg[query_coln], left_on = 'Gene Name', right_on  = 'Worm', how = 'inner')
df_merge.to_csv(os.path.join(fpath_out, fname_worm.replace('.csv','_kf_inner.csv')), index=False)


####### Section 2: Comparison with the Walther 2015 worm paper
## part 1. list all identified aggregate from Walther 2015 Wild-Type
fname_worm = 'Walther2015_TS1D.csv'
fpath_out = 'ProteomeComp/Walther2015_worm'
worm_agg = pd.read_csv(os.path.join(fpath_out, fname_worm))
search_coln = ['log abs.LFQ AgWTD01.1 Corrected','log abs.LFQ AgWTD01.2 Corrected','log abs.LFQ AgWTD01.3 Corrected','log abs.LFQ AgWTD01.4 Corrected',\
'log abs.LFQ AgWTD06 Corrected',\
'log abs.LFQ AgWTD12.1 Corrected','log abs.LFQ AgWTD12.2 Corrected','log abs.LFQ AgWTD12.3 Corrected','log abs.LFQ AgWTD12.4 Corrected',\
'log abs.LFQ AgWTD17 Corrected']
worm_agg.dropna(subset =search_coln, how = 'all', inplace=True)
worm_gene_agg = worm_agg.dropna(subset =['Gene Name'])
worm_gene_agg_tot = worm_gene_agg.shape[0]
print ('total # of proteins identified as aggregate in WT C. elegans in Table S1D is %s'%(worm_agg.shape[0]))
print ('total # of aggregates that can be matched from in Table S1D is %s'%(worm_gene_agg_tot))
# look for overlap with our killifish data
query_coln = ['Tissue', 'Sample.Type', 'N. furzeri Protein Id', 'Human', 'Worm','log2_Young_zscore','log2_Old_zscore','OvY_logFC', 'OvY_pval']
df_merge = worm_gene_agg.merge(kf_agg[query_coln], left_on = 'Gene Name', right_on  = 'Worm', how = 'inner')
df_merge_nodup = df_merge.drop_duplicates(subset =['Gene Name'])
worm_kf_agg_overlap = df_merge_nodup.shape[0]
kf_agg_tot = kf_agg.drop_duplicates(subset=['N. furzeri Protein Id']).shape[0]
print ('total # of aggregates overlap between C. elegan and N. furzeri aggregate is %s'%(worm_kf_agg_overlap))
for tissue in df_merge['Tissue'].unique():
    print ('total # of aggregates in %s is %s, overlap between C. elegan and N. furzeri aggregate is %s'\
           %(tissue, kf_agg[kf_agg['Tissue']==tissue].shape[0], df_merge[df_merge['Tissue']==tissue].shape[0]))
# output files
df_merge.to_csv(os.path.join(fpath_out, fname_worm.replace('.csv','_kf_inner.csv')), index=False)
venn([worm_gene_agg_tot - worm_kf_agg_overlap, kf_agg_tot - worm_kf_agg_overlap, worm_kf_agg_overlap],\
    ['C. elegans', 'N. furzeri'], 'C.elegans_N.furzeri_agg_overlap.pdf', fpath_out)
### part 2. list all identified proteins from Walther 2015 Wild-Type
fname_worm_tl = 'Walther2015_TS1B.csv'
fpath_out = 'ProteomeComp/Walther2015_worm'
worm_tl = pd.read_csv(os.path.join(fpath_out, fname_worm_tl))
search_coln = ['log abs.LFQ WTday06','log abs.LFQ WTday12','log abs.LFQ WTday17','log abs.LFQ WTday22']
worm_tl.dropna(subset =search_coln, how = 'all', inplace=True)
worm_gene_tl = worm_tl.dropna(subset=['Gene Name'])
worm_gene_tl_tot = worm_gene_tl.shape[0]
print ('total # of proteins identified as total protein in WT C. elegans in Table S1B is %s'%(worm_tl.shape[0]))
print ('total # of proteins that can be matched from in Table S1B is %s'%(worm_gene_tl_tot))
# look for overlap with our killifish data
query_coln = ['Tissue', 'Sample.Type', 'N. furzeri Protein Id', 'Human', 'Worm','log2_Young_zscore','log2_Old_zscore','OvY_logFC', 'OvY_pval']
df_merge_tl = worm_gene_tl.merge(kf_tl[query_coln], left_on = 'Gene Name', right_on  = 'Worm', how = 'inner')
df_merge_tl_nodup = df_merge_tl.drop_duplicates(subset =['Gene Name'])
worm_kf_tl_overlap = df_merge_tl_nodup.shape[0]
kf_tl_tot = kf_tl.drop_duplicates(subset=['N. furzeri Protein Id']).shape[0]
print ('total # of overlap between C. elegan and N. furzeri total protein is %s'%(worm_kf_tl_overlap))
for tissue in df_merge_tl['Tissue'].unique():
    print ('total # of proteins in %s is %s, overlap between C. elegan and N. furzeri proteins is %s'\
        %(tissue, kf_tl[kf_tl['Tissue']==tissue].shape[0], df_merge_tl[df_merge_tl['Tissue']==tissue].shape[0]))
# output files
df_merge_tl.to_csv(os.path.join(fpath_out, fname_worm_tl.replace('.csv','_kf_inner.csv')), index=False)
venn([worm_gene_tl_tot - worm_kf_tl_overlap, kf_tl_tot - worm_kf_tl_overlap, worm_kf_tl_overlap],\
    ['C. elegans', 'N. furzeri'], 'C.elegans_N.furzeri_tl_overlap.pdf', fpath_out)
## part 3. compare TL and AGG overlap together
df_worm_tl_agg = worm_gene_agg.merge(worm_gene_tl, left_on = 'Gene Name', right_on = 'Gene Name', how = 'inner')
df_merge_worm_nodup = df_worm_tl_agg.drop_duplicates(subset =['Gene Name'])
worm_tl_agg_overlap = df_merge_worm_nodup.shape[0]
df_merge_tl.to_csv(os.path.join(fpath_out, fname_worm_tl.replace('.csv','_kf_inner.csv')), index=False)
venn([worm_gene_tl_tot - worm_tl_agg_overlap, kf_tl_tot - worm_tl_agg_overlap, worm_tl_agg_overlap],\
    ['C. elegans TL', 'C. elegans AGG'], 'C.elegans_tl_agg_overlap.pdf', fpath_out)


####### Section 3: list all identified aggregate from Ori 2020 Wild-Type
query_coln = ['Tissue', 'Sample.Type', 'N. furzeri Protein Id', 'Human', 'Mouse','log2_Young_zscore','log2_Old_zscore','OvY_logFC', 'OvY_pval']
kf_brain_agg = kf_agg[query_coln][kf_agg['Tissue']=='Brain']
kf_brain_tl = kf_tl[query_coln][kf_tl['Tissue']=='Brain']
#### part 1. comparison against mice aggregate data
fname_orimice = 'EV7_tab1.csv'
fpath_out = 'ProteomeComp/Ori2020_killifish'
# read the Ori mice aggregate data
mice_data = pd.read_csv(os.path.join(fpath_out, fname_orimice))
# expand the gene.name column into rows to make sure i map everything
mice_expand = pd.DataFrame(mice_data['gene.name'].str.split(';').tolist(), index=mice_data.ID).stack()
mice_expand = mice_expand.reset_index([0, 'ID'])
mice_expand.columns = ['ID', 'gene.name']
mice_expand = mice_expand.merge(mice_data, left_on='ID', right_on='ID', how='left',suffixes=['', '_old'])
mice_tot = mice_data.shape[0]
print ('total # of proteins identified as aggregate in mice brain in Table EV7.tab1 is %s'%(mice_tot))
# look for overlap with our data
df_merge = mice_expand.merge(kf_brain_agg, left_on = 'gene.name', right_on  = 'Mouse', how = 'inner')
df_merge_nodup = df_merge.drop_duplicates(subset =['gene.name'])
mice_kf_agg_overlap = df_merge_nodup.shape[0]
kf_agg_tot = kf_brain_agg.shape[0]
print ('total # of aggregates overlap between mice and N. furzeri brain is %s'%(mice_kf_agg_overlap))
# output files
df_merge.to_csv(os.path.join(fpath_out, fname_orimice.replace('.csv','_kfbrainagg_inner.csv')), index=False)
venn([mice_tot - mice_kf_agg_overlap, kf_agg_tot - mice_kf_agg_overlap, mice_kf_agg_overlap],\
    ['Mice brain ori', 'N. furzeri brain'], 'Mice_N.furzeri_brain_agg_overlap.pdf', fpath_out)
#### part 2. comparison against killifish aggregate data
fname_orikf = 'EV7_tab2.csv'
fpath_out = 'ProteomeComp/Ori2020_killifish'
# read the Ori aggregate data
ori_data = pd.read_csv(os.path.join(fpath_out, fname_orikf))
ori_tot = ori_data.shape[0]
print ('total # of proteins identified as aggregate in Ori et al. killifish brain in Table EV7.tab2 is %s'%(ori_tot))
# look for overlap with our data
df_merge = ori_data.merge(kf_brain_agg, left_on = 'Human_gene_name', right_on  = 'Human', how = 'inner')
df_merge_nodup  = df_merge.drop_duplicates(subset='Human')
kf_agg_overlap = df_merge_nodup.shape[0]
kf_agg_tot = kf_brain_agg.shape[0]
print ('total # of aggregates overlap between Ori and our (N. furzeri brain) data is %s'%(kf_agg_overlap))
# output files
df_merge.to_csv(os.path.join(fpath_out, fname_orikf.replace('.csv','_kfbrainagg_inner.csv')), index=False)
venn([ori_tot - kf_agg_overlap, kf_agg_tot - kf_agg_overlap, kf_agg_overlap],\
    ['Killifish brain ori', 'N. furzeri brain'], 'Ori_Our_brain_agg_overlap.pdf', fpath_out)
#### part 3. comparison against killifish tissue lysate data
### i need actual map to convert their protein ID to NCBI ID, skip this part! ####
fname_orikf = 'EV2_tab2.csv'
fpath_out = 'ProteomeComp/Ori2020_killifish'
# read the Ori aggregate data
ori_data = pd.read_csv(os.path.join(fpath_out, fname_orikf))
ori_tot = ori_data.shape[0]
print ('total # of proteins identified in total proteome in Ori et al. killifish brain in Table EV2.tab2 is %s'%(ori_tot))
# look for overlap with our data
df_merge = ori_data.merge(kf_brain_tl, left_on = 'human.gene.name', right_on  = 'Human', how = 'inner')
df_merge_nodup  = df_merge.drop_duplicates(subset='Human')
kf_tl_overlap = df_merge_nodup.shape[0]
kf_tl_tot = kf_brain_agg.shape[0]
print ('total # of proteins overlap between Ori and our (N. furzeri brain) data is %s'%(kf_agg_overlap))
# output files
df_merge.to_csv(os.path.join(fpath_out, fname_orikf.replace('.csv','_kfbraintl_inner.csv')), index=False)
venn([ori_tot - kf_agg_overlap, kf_agg_tot - kf_agg_overlap, kf_agg_overlap],\
    ['Killifish brain ori', 'N. furzeri brain'], 'Ori_Our_brain_tl_overlap.pdf', fpath_out)


####### Section 4: list all identified aggregate from David 2010 Table S1
query_coln = ['Tissue', 'Sample.Type', 'N. furzeri Protein Id', 'Human', 'Worm','log2_Young_zscore','log2_Old_zscore','OvY_logFC', 'OvY_pval']
fname_worm = 'David_TableS1.csv'
fpath_out = 'ProteomeComp/David2010_kenyon'
# import the worm data
worm_data = pd.read_csv(os.path.join(fpath_out, fname_worm))
worm_data.dropna(subset =['Gene name'], inplace=True)
worm_tot = worm_data.shape[0]
print ('total # of proteins identified as aggregate in worm from David et al. 2009 in Table S1 is %s'%(worm_tot))
df_merge = worm_data.merge(kf_agg[query_coln], left_on = 'Gene name', right_on  = 'Worm', how = 'inner')
df_merge_nodup = df_merge.drop_duplicates(subset =['Gene name'])
worm_kf_agg_overlap = df_merge_nodup.shape[0]
kf_agg_tot = kf_agg.drop_duplicates(subset=['N. furzeri Protein Id']).shape[0]
print ('total # of aggregates overlap between worm and N. furzeri is %s'%(worm_kf_agg_overlap))
# age associated aggregate
df_merge['Worm.Killifish.UpWithAge'] = df_merge.apply(David_age_agg, axis=1)
worm_kf_age_agg = df_merge['Gene name'][df_merge['Worm.Killifish.UpWithAge']=='YES'].drop_duplicates()
print ('total # of age-associated aggregates (increase with age) that overlap between worm and N. furzeri is %s'\
       %(worm_kf_age_agg.shape[0]))
# output files
df_merge.to_csv(os.path.join(fpath_out, fname_worm.replace('.csv','_kfagg_inner.csv')), index=False)
venn([worm_tot - worm_kf_agg_overlap, kf_agg_tot - worm_kf_agg_overlap, worm_kf_agg_overlap],\
    ['C. elegans David 2010', 'N. furzeri'], 'Worm_N.furzeri_agg_overlap.pdf', fpath_out)


###### Section 5: draw a venn3 diagram to compare the age-associated aggregate among worm and killifish dataset
A = 'David 2010'
fpath_A = 'ProteomeComp/David2010_kenyon'
fname_A = 'David_TableS1.csv'
df_A = pd.read_csv(os.path.join(fpath_A, fname_A))
A_idx = np.where(df_A['TwoExp>1.5']=="YES")
A_worm = df_A['Gene name'].loc[A_idx].dropna()

B = 'Walther 2015'
fpath_B = 'ProteomeComp/Walther2015_worm'
fname_B = 'Walther2015_TS1D.csv'
df_B = pd.read_csv(os.path.join(fpath_B, fname_B))
B_idx = np.where(df_B['log abs.LFQ AgWTD17 Corrected']> df_B['log abs.LFQ AgWTD06 Corrected'])
B_worm = df_B['Gene Name'].loc[B_idx].dropna()

C = 'N furzeri 2020'
C_idx = np.where((kf_agg['OvY_logFC']>=0) & (kf_agg['OvY_pval']<=0.05))
df_C =  kf_agg['N. furzeri Protein Id'].iloc[C_idx].drop_duplicates()
C_worm = kf_agg['Worm'].iloc[C_idx].dropna()

venn([set(A_worm), set(B_worm), set(C_worm)], [A,B,C], 'Worm_ageagg.pdf','ProteomeComp')
