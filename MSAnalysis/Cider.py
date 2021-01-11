#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 14:46:14 2018

@author: yiwenchen

make sure localCIDER package is installed, see https://github.com/Pappulab/localCIDER
"""

# this script is used to run Cider developed by the Pappu lab with the loalCider package
import localcider
from localcider.sequenceParameters import SequenceParameters
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import datetime
import multiprocessing as mp


#define a function where you specify sequence and output all properties in a dictionary
def cider_seq(seq, seq_name, seq_XP_id):
    seq_dict = {'Seq': seq, 'Name': seq_name, 'Protein Id': seq_XP_id}
    #import an amino acid sequence to create the SequenceParameters object
    SeqOb=SequenceParameters(seq)
    #output integer equal to the sequence length
    seq_dict['length'] = SeqOb.get_length()
    #output a protein's mean hydropathy (a value between 0 and 9, where 0 is the least hydrophobic and 9 is the most). This is simply a re-based Kyte-Doolittle scale, which runs from 0 to 9 instead of from -4.5 to 4.5, as in the original paper.
    seq_dict['mean_hydropathy'] = SeqOb.get_mean_hydropathy()
    #Get a protein's mean hydropathy as defined by Uversky, using a normalized Kyte-Doolittle hydropathy table. Here, values range between 0 and 1, with 0 being the least hydrophobic and 1 the most hydrophobic.
    seq_dict['uversky_hydropathy'] = SeqOb.get_uversky_hydropathy()
    #Get a protein's fraction of residues which are considered '[D]isorder promoting'
    seq_dict['FracDisoPromoting'] = SeqOb.get_fraction_disorder_promoting()
    #Return a dictionary with the fractions of each amino acid in your sequence
    seq_AAFrac = SeqOb.get_amino_acid_fractions()
    for i in seq_AAFrac.keys():
        seq_dict['Frac_%s'%i] = seq_AAFrac[i]
    #Get the sequence charge decoration (SCD) value associated with a sequence.
    #SeqOb.get_SCD()
    #Get the kappa value associated with a sequence.
    seq_dict['kappa'] = SeqOb.get_kappa()
    #Get the Omega value associated with a sequence. Omega describes the patterning between charged/proline residues and all other residues.  
    seq_dict['Omega'] = SeqOb.get_Omega()
    #Get the 2-alphabet sequence used for calculating the Omega parameter as defined in Martin et al. R/K/D/E/P are represented as X and all other residues are O.
    seq_dict['Omega_seq'] = SeqOb.get_Omega_sequence()
    #Get the user defined patterning paremeter, where residues are grouped into either two groups (grp1 and ALL other residues) or three groups (grp1, grp2, and ALL other residues).
    seq_dict['kappa_X'] = SeqOb.get_kappa_X(['E','D'], grp2 = ['K','D'])
    #Get the maximum delta value for a sequence of this composition. Note kappa is delta/deltaMax.
    seq_dict['deltaMax'] = SeqOb.get_deltaMax()
    #Get the delta value for this specific sequence. Note kappa is delta/deltaMax.
    seq_dict['delta'] = SeqOb.get_delta()
    #Get the number of positive residues in the sequence
    seq_dict['countPos'] = SeqOb.get_countPos()
    #Get the number of negative residues in the sequence
    seq_dict['countNeg'] = SeqOb.get_countNeg()
    #Get the number of neutral residues in the sequence
    seq_dict['countNeut'] = SeqOb.get_countNeut()
    #Get the fraction of positively charged residues in the sequence
    seq_dict['FracPos'] = SeqOb.get_fraction_positive()
    #Get the fraction of negatively residues in the sequence
    seq_dict['FracNeg'] = SeqOb.get_fraction_negative()
    #Get the fraction of charged residues in the sequence.
    seq_dict['FracCharged'] = SeqOb.get_FCR() #can change individual AA's pH contribution
    #Get the fraction of expanding residues in the sequence. We define 'expanding' residues as D/E/R/K/P. This will be the same as the FCR+fraction of proline
    seq_dict['FracExpanding'] = SeqOb.get_fraction_expanding()
    #Get the net charge per residue of the sequence.
    seq_dict['NCPR'] = SeqOb.get_NCPR()
    #Get the absolute magnitude of the mean net charge
    seq_dict['MeanNetCharge'] = SeqOb.get_mean_net_charge()
    #Get the absolute magnitude of the mean net charge
    seq_dict['pI'] = SeqOb.get_isoelectric_point()
    #Returns the molecular weight of the protein in Da
    seq_dict['MW'] = SeqOb.get_molecular_weight()
    #Returns the IDP diagram of states [REF 1] region based on the FCR/NCPR rules. 
    seq_dict['phasePlotRegion'] = SeqOb.get_phasePlotRegion()
    #Function which recomputes kappa after complete phosphorylation based on the currently defined phosphosites.
    seq_dict['kappa_afterPhos'] = SeqOb.get_kappa_after_phosphorylation()
    #Function used to conduct position-specific sequence analysis
    seq_dict['FCR_blobLen5'] = ','.join(SeqOb.get_linear_FCR(blobLen=5)[1].astype(str))
    seq_dict['NCPR_blobLen5'] = ','.join(SeqOb.get_linear_NCPR(blobLen=5)[1].astype(str))
    seq_dict['sigma_blobLen5'] = ','.join(SeqOb.get_linear_sigma(blobLen=5)[1].astype(str))
    seq_dict['hydropathy_blobLen5'] = ','.join(SeqOb.get_linear_hydropathy(blobLen=5)[1].astype(str))
    return seq_dict


#define a function to generate the phase plot using scatter plot
def phaseplot_scatter(x, y, xlabel, ylabel, title, figtitle, alpha = 0.1):
    fig, ax = plt.subplots()
    ax.scatter(x, y, alpha = alpha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_title(title)
    fig.savefig(figtitle)
    
  
#define a function to generate the phase plot using 2d histogram plot
def phaseplot_2dhist(x, y, xlabel, ylabel, title, figtitle, cmap = plt.cm.jet, bins = 50):
    fig, ax = plt.subplots()
    ax.hist2d(x, y, bins=(bins, bins), cmap=cmap)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_title(title)
    fig.savefig(figtitle)
 
    
#define a function to generate the phase plot using seaborn jointplot (scatter/kde) & histogram
def phaseplot_joint(x, y, df, xlabel, ylabel, title, figtitle, figtype = 'scatter', margin = False,alpha=0.05):
    g = sns.JointGrid(x=x, y=y, data=df, xlim=(0,1), ylim=(0,1),space=0)
#    g.plot_joint(plt.scatter, color='k', alpha = 0.1)
#    g.plot_marginals(sns.distplot, kde=False, color=".5")
    g.ax_joint.fill([0,0.35,0], [1,0.675 ,0.35], color = 'r', label= 'Negatively charged strong polyelectrolytes:\nSwollen coils')
    g.ax_joint.fill([0.35,0.675,1], [0,0.35 ,0], color = 'b', label='Positively charged strong polyelectrolytes:\nSwollen coils')
    g.ax_joint.fill([0,0.35,0.675, 0.35], [0.35,0.675,0.35,0], color = 'g', label='Strong polyampholytes:\nCoils, hairpins, & chimeras')
    g.ax_joint.fill([0,0,0.35,0.25], [0.25,0.35,0,0], color = '#66c2a4', label='Janus sequences:\nCollapsed or expanded - context dependent')
    g.ax_joint.fill([0,0,0.25], [0, 0.25,0], color = '#d9f0a3', label='Weak polyampholytes & polyelectrolytes:\nGlobules & tadpoles')
    g.ax_joint.set_xlabel(xlabel)
    g.ax_joint.set_ylabel(ylabel)        
    if figtype == 'scatter':
        g.plot_joint(plt.scatter, color='k', alpha = alpha, zorder = 10)
    elif figtype == 'kde':
        g.plot_joint(sns.kdeplot, color='k', zorder = 10)
    if margin == True:
        g.plot_marginals(sns.distplot, kde=True, color='k')
    g.fig.suptitle(title)
    g.ax_joint.legend()
    g.savefig(figtitle)


# updated on 2020.04.02     
def run_cider(df, seq_coln, name_coln, id_coln, fpath_out, fname):
    df_min = df[[seq_coln, name_coln, id_coln]].dropna(subset = [seq_coln, id_coln], how='any')
    all_ciderOb = df_min.apply(lambda x: cider_seq(x[seq_coln], x[name_coln], x[id_coln]), axis=1)
    df_cider = pd.DataFrame(list(all_ciderOb))
    final_coln = ['Name', 'Protein Id', 'Seq']#, 'MW', 'pI','length', 'mean_hydropathy', 'uversky_hydropathy']
    for i in df_cider.columns:
        if i not in final_coln: final_coln.append(i) 
    df_cider[final_coln].to_csv(os.path.join(fpath_out, fname.replace('.csv','_Cider.csv')),index = False)
    print (datetime.datetime.now())
    return 'Done'########### SECTION 11: This is used to run cider function in parallel  ##########


# ############ SECTION 1: THIS IS USED TO FEED THE WHOLE KILLFISH PROTEOME TO CIDER##########
# #################### Use Parallel computing to re-run killifish cider ####################
# ## updated on 2020.04.13, note that code will take couple hours to finish!!!      
# ########## Split files and do parallel computing ######
# begin_time = datetime.datetime.now()
# fpath_proteome = 'Properties'
# fpath_out = 'Properties/cider'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# n_sub = 100
# df_proteome = pd.read_csv(os.path.join(fpath_proteome, 'Killifish_DISOPRED_PLAAC.csv'))
# ###### split a single dataframe into multiple ones ######
# df_array = np.array_split(df_proteome, n_sub)
# for i, df_i in enumerate(df_array):
#     df_i.to_csv(os.path.join(fpath_out, 'S%s.csv')%i, index=False)
# ########## parallel computing ####
# fname_prefix = 'S'
# custom_range = range(49, n_sub)
# fnames = ['%s%s.csv'%(fname_prefix, i) for i in custom_range]
# f_list = [(pd.read_csv(os.path.join(fpath_out, f)), f) for f in fnames]
# arg_list = [(df_i, 'Seq', 'Name', 'Protein Id', fpath_out, fname_i) for df_i, fname_i in f_list]
# pool = mp.Pool(mp.cpu_count())
# results = [pool.apply_async(run_cider, args = arg_i) for arg_i in arg_list]
# pool.close()
# print(datetime.datetime.now() - begin_time)
# ###### end of parallel computing ####
# # stich the result files together
# for i in range(50):
#     fname_i = '%s%s_Cider.csv'%(fname_prefix, i)
#     df_i = pd.read_csv(os.path.join(fpath_out, fname_i))
#     if i == 0:
#         df = df_i
#     else:
#         df = df.append(df_i, ignore_index=True)
# df.to_csv(os.path.join(fpath_out, 'Killifish_Cider.csv'),index=False)
# print(datetime.datetime.now() - begin_time)
# ###################################### End of Section 11 ####################################


############ SECTION 2: THIS IS USED TO FEED THE WHOLE KILLFISH PROTEOME TO CIDER##########
# use this to run a selectiv few seq for code check
#import Killifish DISO_PLAAC file    
fpath_proteome = 'Properties'
fpath_out = 'Properties/cider'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
df_proteome = pd.read_csv(os.path.join(fpath_proteome, 'Killifish_DISOPRED_PLAAC.csv'))
#feed in sequences one at a time and propogate each proteins's cider output to a dataframe
all_ciderOb = []
excluded_aa = ['B', 'U', 'O', 'X', 'J', 'Z']
### specify the rows you want to run by index number
df_sample = df_proteome.iloc[1:5,:]
for index, row in df_sample[['Seq','Name','Protein Id']].iterrows():
    seq_name = row['Name']
    seq_XP_id = row['Protein Id']
    seq = row['Seq']
    if len(seq) !=0:
        if len([i for i in excluded_aa if i in seq])==0:
            row_dict = cider_seq(seq, seq_name, seq_XP_id)
            all_ciderOb.append(row_dict)
        else:
            print (seq_XP_id)
    else:
        print (seq_XP_id)
df_proteome_cider = pd.DataFrame(all_ciderOb)
final_coln = ['Human', 'N. furzeri Protein Id', 'Seq', 'MW', 'pI','length', 'mean_hydropathy', 'uversky_hydropathy']
for i in df_proteome_cider.columns:
    if i not in final_coln: final_coln.append(i) 
df_proteome_cider[final_coln].to_csv(os.path.join(fpath_out, 'Killifish_Cider_example.csv'),index = False)
##################################### End of Section 1 ####################################


# ############ SECTION 3: Run a few seqs of interest ##########
# fpath_proteome = 'Properties'
# fpath_out = 'Properties/cider'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# df_proteome = pd.read_csv(os.path.join(fpath_proteome, 'HumanPKM.csv'))
# all_ciderOb = []
# excluded_aa = ['B', 'U', 'O', 'X', 'J', 'Z']
# ### specify the rows you want to run by index number
# df_sample = df_proteome
# for index, row in df_sample[['Seq','Name','Protein Id']].iterrows():
#     seq_name = row['Name']
#     seq_XP_id = row['Protein Id']
#     seq = row['Seq']
#     if len(seq) !=0:
#         if len([i for i in excluded_aa if i in seq])==0:
#             row_dict = cider_seq(seq, seq_name, seq_XP_id)
#             all_ciderOb.append(row_dict)
#         else:
#             print (seq_XP_id)
#     else:
#         print (seq_XP_id)
# df_proteome_cider = pd.DataFrame(all_ciderOb)
# final_coln = ['Name', 'Protein Id', 'Seq', 'MW', 'pI','length', 'mean_hydropathy', 'uversky_hydropathy']
# for i in df_proteome_cider.columns:
#     if i not in final_coln: final_coln.append(i) 
# df_proteome_cider[final_coln].to_csv(os.path.join(fpath_out, 'HumanPKM_Cider.csv'),index = False)
# ##################################### End of Section 2 ####################################


############# SECTION 4: This is used to save and load object saved to file #############
#def save_obj(obj, name ):
#    with open('obj/'+ name + '.pkl', 'wb') as f:
#        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
#
#def load_obj(name ):
#    with open('obj/' + name + '.pkl', 'rb') as f:
#        return pickle.load(f)
#################################### End of Section 4 ####################################
    

########### SECTION 5: This is used to extract all parameters for SysSeedingProteins ##########
#seeding_fname = 'SysSeedingProteinsRaw.csv'
#seeding_fpath = 'Correlomics/SysSeeding'
#database_fname = 'Killifish_DISOPRED_PLAAC.csv'
#database_fpath = 'Correlomics'
#cider_fname = 'Killifish_Cider.csv'
#cider_fpath = 'Correlomics'
#df_seeding = pd.read_csv(os.path.join(seeding_fpath, seeding_fname))
#df_database = pd.read_csv(os.path.join(database_fpath, database_fname))
#df_cider = pd.read_csv(os.path.join(cider_fpath, cider_fname))
#df_all = pd.merge(df_seeding, df_database, how = 'left', left_on = ['Protein Id','seq'], right_on = ['Protein Id','Seq'])
#df_all = pd.merge(df_all, df_cider, how = 'left', left_on = [], right_on = ['Protein Id','Seq'])
#df_all.to_csv(os.path.join(seeding_fpath, seeding_fname.replace('.csv', 'DISOLLRCider.csv')))
#
############ SECTION 6: This is used to extract all parameters for SysSeedingProteins ##########
#seeding_fname = 'SysSeedingProteinsRaw.csv'
#seeding_fpath = 'Correlomics/SysSeeding'
#database_fname = 'Killifish_DISOPRED_PLAAC.csv'
#database_fpath = 'Correlomics'
#df_seeding = pd.read_csv(os.path.join(seeding_fpath, seeding_fname))
#df_database = pd.read_csv(os.path.join(database_fpath, database_fname))
#df_all = pd.merge(df_seeding, df_database, how = 'left', left_on = ['N. furzeri Protein Id'], right_on = ['Protein Id'])
###import Killifish DISO_PLAAC file    
##in_path = 'Correlomics/SysSeeding'
##in_fname = 'SysSeedingProteinsRaw.csv'
##df_proteome = pd.read_csv(os.path.join(in_path, in_fname))
##feed in sequences one at a time and propogate each proteins's cider output to a dataframe
#all_ciderOb = []
#excluded_aa = ['B', 'U', 'O', 'X', 'J', 'Z']
#for index, row in df_all[['Seq','Human','N. furzeri Protein Id']].iterrows():
#    seq_name = row['Human']
#    seq_XP_id = row['N. furzeri Protein Id']
#    seq = row['Seq']
#    if len(seq) !=0:
#        if len([i for i in excluded_aa if i in seq])==0:
#            row_dict = cider_seq(seq, seq_name, seq_XP_id)
#            all_ciderOb.append(row_dict)
#        else:
#            print seq_XP_id
#    else:
#        print seq_XP_id
#df_proteome_cider = pd.DataFrame(all_ciderOb)
#final_coln = ['Human', 'N. furzeri Protein Id', 'Seq', 'MW', 'pI','length', 'mean_hydropathy', 'uversky_hydropathy']
#for i in df_proteome_cider.columns:
#    if i not in final_coln: final_coln.append(i) 
#df_proteome_cider[final_coln].to_csv(os.path.join(seeding_fpath, seeding_fname.replace('.csv', '_Cider.csv')),index = False)
#df_combine = pd.merge(df_all, df_proteome_cider, how = 'left', left_on = ['N. furzeri Protein Id', 'Human', 'Seq'], right_on = ['N. furzeri Protein Id', 'Human', 'Seq'])
#df_combine.to_csv(os.path.join(seeding_fpath, seeding_fname.replace('.csv', '_DISOLLRCider.csv')),index = False)