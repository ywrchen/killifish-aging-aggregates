#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:31:30 2019

# this file is used to plot geneset heatmap in bubble plot style

@author: yiwenchen
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
import matplotlib as mpl


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 20})
sns.set(style="white")


# function used to generate bubble plot for table with two variable of interest
# the input df needs to be flat where each row presents a single pair of var1 and var2 for a given protein
def bubble(df, var_x, var_y, var_value, var_size, fname, path, colormap = 'default', fig_w = '', fig_l = '',\
           var_x_order = 'None', var_y_order = 'None', order = 'AsIs', var_size_cutoff = 0.05, xaxis_rot = 0, yaxis_rot = 0):
    # var_x will be presented on x axis while var_y will be presented on y axis
    # example: var_x can be tissue, var_y can be Human, var_value can be YvO_logFC, var_size can be YvO_pval
    if var_x_order == 'None':
        x_unique, X = str_to_num(df[var_x].tolist(), order = 'AsIs')
        x_dict = dict_twolist(X, df[var_x])
    else:
        x_unique = var_x_order.keys()
        X = [var_x_order[i] for i in df[var_x].tolist()]
        x_dict = dict_twolist([var_x_order[i]for i in x_unique], x_unique) #dict_twolist(X, df[var_x])
    if var_y_order == 'None':
        y_unique, Y = str_to_num(df[var_y].tolist(), order = order)
    elif var_y_order == 'classification':
        df_ylabels = df[['classification',var_y]].sort_values(by=['classification',var_y])
        y_unique = df_ylabels[var_y].drop_duplicates().tolist()
        if order == 'AsIsR':
            Y = [(len(y_unique)-y_unique.index(i)-1) for i in df[var_y].tolist()]
        else:
            Y = [y_unique.index(i) for i in df[var_y].tolist()]
    else:
        y_unique = var_y_order.keys()
        Y = [var_y_order[i] for i in df[var_y].tolist()]
    y_dict = dict_twolist(Y, df[var_y])
    if fig_w == '': fig_w = len(x_unique)
    if fig_l == '': fig_l = len(y_unique)*0.65
#    if fig_l <6: fig_l = 6
    # if you want the y be in ascending order as the human term go from a->z in asending ordering as well
    Z = df[var_value].tolist()
    S = -np.log10(df[var_size])
    df = pd.DataFrame(list(zip(X,Y,S,Z)), columns = ['X','Y','S','Z'])
    df_sig = df[df['S'] > -np.log10(var_size_cutoff)]
    df_nosig = df[df['S'] <= -np.log10(var_size_cutoff)]
    if colormap == 'default':
        colormap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])   
        cmap_vmin = -0.5; cmap_vmax = 0.5
        cnorm = mpl.colors.Normalize(cmap_vmin, cmap_vmax)
    fig, ax = plt.subplots(1,1, figsize=(fig_w, fig_l))
    size_scale = 400
    ax.scatter(df_nosig['X'], df_nosig['Y'], s = df_nosig['S']*size_scale, c = df_nosig['Z'], cmap = colormap, \
               norm = cnorm, marker = '.', edgecolor = 'grey', linewidth = 0.75, label='_nolegend_') #marker = 'X'
    sc = ax.scatter(df_sig['X'], df_sig['Y'], s = df_sig['S']*size_scale, c = df_sig['Z'], cmap = colormap, \
                    norm = cnorm, marker = 'o', edgecolor = 'grey', linewidth = 0.75, label='_nolegend_')
    ax.set_title(fname.replace('_bubble.pdf', ''))
    ax.set(xticks = range(len(x_unique)), xticklabels = [x_dict[x_i] for x_i in range(len(x_unique))],\
           yticks = range(len(y_unique)), yticklabels = [y_dict[y_i] for y_i in range(len(y_unique))])
    # write a legend with the desired size and label
    pval_labels = [1, 0.5, 0.1, 0.05, 0.005]
    pval_sizes = [-np.log10(i) if i!=1 else 0.004364 for i in pval_labels]
    for pval_size, pval_label in zip(np.array(pval_sizes)*size_scale, pval_labels):
        if pval_label >var_size_cutoff:
            ax.scatter([], [], marker = '.', c='k', alpha = 0.5, s=pval_size, label = pval_label)
        else:
            ax.scatter([], [], marker = 'o', c='k', alpha = 0.5, s=pval_size, label = pval_label)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1.13, title="pval", markerscale = 1,
               loc = 'lower left', fontsize = 'large', bbox_to_anchor=(1, 0.1), borderaxespad=0)
    ax.tick_params(axis='x', labelrotation=xaxis_rot)
    ax.tick_params(axis='y', labelrotation=yaxis_rot)
#    # place the fold change colorbar
    cbaxes = fig.add_axes([0.93, 0.5, 0.03, 0.3])
    plt.colorbar(sc, shrink= 0.5, aspect = 10, cax = cbaxes)
    #    plt.tight_layout()
    fig.savefig(os.path.join(path, fname), bbox_inches = 'tight')
    return None


# function used to generate bubble plot for table with two variable of interest
# the input df is in the wide form, basically need to convert the matrix into scatter plot and preserve the row and column index
# the x axis labbel is the the column name, the y axis label is just the row index from the input dataframe
# the color of each point is specified in df_value, and the sizes are specified in df_size; both df have to have identical dimensions and match 1:1
def bubble_wide(df_value, df_size, fname_fig, fpath_fig, colormap = 'PiYG', fig_w = '', fig_l = '', size_scale = 50,\
                size_cutoff = 0.05, xaxis_rot = 0, yaxis_rot = 0, y_order = 'as_matrix', cmap_vmin = -1, cmap_vmax = 1): 
    ylim, xlim = df_value.shape
    X = range(xlim)*ylim # used to track the row index
    if y_order == 'as_matrix':
        Y = np.repeat(range(ylim-1,-1,-1),xlim) # used to track the column index
        ylabels = df_value.index[::-1]
    elif y_order == 'reverse':
        Y = np.repeat(range(ylim),xlim) # used to track the column index
        ylabels = df_value.index        
    S = -np.log10(df_size.values).flatten() # used to track the size of the scatter points
    Z = df_value.values.flatten() # used to track the value/color of the scatter points
    if fig_w == '': fig_w = xlim
    if fig_l == '': fig_l = ylim*0.65
#    if fig_l <6: fig_l = 6
    # if you want the y be in ascending order as the human term go from a->z in asending ordering as well
    df = pd.DataFrame(list(zip(X,Y,S,Z)), columns = ['X','Y','S','Z'])
    df_sig = df[df['S'] > -np.log10(size_cutoff)]
    df_nosig = df[df['S'] <= -np.log10(size_cutoff)]
    cmap_vmin = cmap_vmin; cmap_vmax = cmap_vmax
    if colormap == 'default':
        colormap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])   
        cnorm = mpl.colors.Normalize(cmap_vmin, cmap_vmax)
    else:
        colormap = mpl.cm.get_cmap(colormap)
        cnorm = mpl.colors.Normalize(cmap_vmin, cmap_vmax)
    sns.set(style="white")
    fig, ax = plt.subplots(1,1, figsize=(fig_w, fig_l))
    size_scale = size_scale
    ax.scatter(df_nosig['X'], df_nosig['Y'], s = df_nosig['S']*size_scale, c = df_nosig['Z'], cmap = colormap, \
               norm = cnorm, marker = '.', edgecolor = 'grey', linewidth = 0.75, label='_nolegend_') #marker = 'X'
    sc = ax.scatter(df_sig['X'], df_sig['Y'], s = df_sig['S']*size_scale, c = df_sig['Z'], cmap = colormap, \
                    norm = cnorm, marker = 'o', edgecolor = 'grey', linewidth = 0.75, label='_nolegend_')
    ax.set_title(fname_fig.replace('_bubble.pdf', ''))
    ax.set(xticks = range(xlim), xticklabels = df_value.columns,\
           yticks = range(ylim), yticklabels = ylabels)
    # write a legend with the desired size and label
    pval_labels = [1, 0.5, 0.1, 0.05, 0.005]
    pval_sizes = [-np.log10(i) if i!=1 else 0.004364 for i in pval_labels]
    for pval_size, pval_label in zip(np.array(pval_sizes)*size_scale, pval_labels):
        if pval_label >size_cutoff:
            ax.scatter([], [], marker = '.', c='k', alpha = 0.5, s=pval_size, label = pval_label)
        else:
            ax.scatter([], [], marker = 'o', c='k', alpha = 0.5, s=pval_size, label = pval_label)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=0.5, title='pval', fontsize=24, markerscale = 1,
               loc = 'lower left', bbox_to_anchor=(1, 0.1), handletextpad =0.05, borderaxespad=0) 
    ax.tick_params(axis='x', labelrotation=xaxis_rot,labelsize=24)
    ax.tick_params(axis='y', labelrotation=yaxis_rot,labelsize=24)
#    # place the fold change colorbar
    cbaxes = fig.add_axes([0.93, 0.5, 0.03, 0.3])
    cb = plt.colorbar(sc, shrink= 0.5, aspect = 10, cax = cbaxes)
    cb.ax.tick_params(labelsize=24)
    #    plt.tight_layout()
    fig.savefig(os.path.join(fpath_fig, fname_fig), bbox_inches = 'tight')
    plt.close(fig)
    return None


# function used to generate bubble plot for table with two variable of interest
# the input df will be in short form where all values are presented in columns already
# this will also filter out fil_index_coln where the row value min is larger than var_size_cutoff (useful for filter by p-value)  
# note you CANNOT specify var_x_order or var_y_order here because of the filetering (KEY MATCHING ERROR)
def bubble_filter(df, var_x, var_y, var_value, var_size, fname, path, fil_index_coln, \
                  var_x_order = 'None', var_y_order = 'None', colormap = 'default', \
                  fig_w = '', fig_l = '', order = 'AsIs', var_size_cutoff = 0.05,\
                  xaxis_rot = 0, yaxis_rot = 0):
    # var_x will be presented on x axis (tisue) while var_y will be presented on y axis (protein)
    # example: var_x can be tissue, var_y can be Human, var_value can be YvO_logFC, var_size can be YvO_pval
    df_pfil = df.groupby([fil_index_coln]).min()
    drop_row = df_pfil.index[df_pfil[var_size] > var_size_cutoff].tolist()
    df = df[~df[fil_index_coln].isin(drop_row)]
    if df.shape[0] == 0:
        return '%s is an empty matrix after filtering' %fname
    if var_x_order == 'None':
        x_unique, X = str_to_num(df[var_x].tolist(), order = 'AsIs')
        x_dict = dict_twolist(X, df[var_x])
    else:
        x_unique = var_x_order.keys()
        X = [var_x_order[i] for i in df[var_x].tolist()]
        x_dict = dict_twolist([var_x_order[i]for i in x_unique], x_unique) #dict_twolist(X, df[var_x])
    if var_y_order == 'None':
        y_unique, Y = str_to_num(df[var_y].tolist(), order = order)
    elif var_y_order == 'classification':
        df_ylabels = df[['classification',var_y]].sort_values(by=['classification',var_y])
        y_unique = df_ylabels[var_y].drop_duplicates().tolist()
        if order == 'AsIsR':
            Y = [(len(y_unique)-y_unique.index(i)-1) for i in df[var_y].tolist()]
        else:
            Y = [y_unique.index(i) for i in df[var_y].tolist()]
    else:
        y_unique = var_y_order.keys()
        Y = [var_y_order[i] for i in df[var_y].tolist()]
    y_dict = dict_twolist(Y, df[var_y])
    if fig_w == '': fig_w = len(x_unique)
    if fig_l == '': fig_l = len(y_unique)*0.65
#    if fig_l <6: fig_l = 6
    # if you want the y be in ascending order as the human term go from a->z in asending ordering as well
#    y_unique, Y = np.unique(df[var_y], return_inverse = True)
    Z = df[var_value].tolist()
    S = -np.log10(df[var_size])
    df_fil = pd.DataFrame(list(zip(X,Y,S,Z)), columns = ['X','Y','S','Z'])
    df_sig = df_fil[df_fil['S'] > -np.log10(var_size_cutoff)]
    df_nosig = df_fil[df_fil['S'] <= -np.log10(var_size_cutoff)]
    if colormap == 'default':
        colormap = mpl.colors.LinearSegmentedColormap.from_list("", ["blue", "white", "yellow"])   
        cmap_vmin = -0.5; cmap_vmax = 0.5
        cnorm = mpl.colors.Normalize(cmap_vmin, cmap_vmax)
    fig, ax = plt.subplots(1,1, figsize=(fig_w, fig_l))
    size_scale = 500
    ax.scatter(df_nosig['X'], df_nosig['Y'], s = df_nosig['S']*size_scale, c = df_nosig['Z'], cmap = colormap, \
               norm = cnorm, marker = '.', edgecolor = 'grey', linewidth = 0.75, label='_nolegend_') #marker = 'X'
    sc = ax.scatter(df_sig['X'], df_sig['Y'], s = df_sig['S']*size_scale, c = df_sig['Z'], cmap = colormap, \
                    norm = cnorm, marker = 'o', edgecolor = 'black', linewidth = 2, label='_nolegend_')
    ax.set_title(fname.replace('_bubble.pdf', ''))
    ax.set(xticks = range(len(x_unique)), xticklabels = [x_dict[x_i] if x_i in x_dict.keys() else '' for x_i in range(len(x_unique))],\
           yticks = range(len(y_unique)), yticklabels = [y_dict[y_i] for y_i in range(len(y_unique))])
    # write a legend with the desired size and label
    pval_labels = [1, 0.5, 0.1, 0.05, 0.005]
    pval_sizes = [-np.log10(i) if i!=1 else 0.004364 for i in pval_labels]
    for pval_size, pval_label in zip(np.array(pval_sizes)*size_scale, pval_labels):
        if pval_label >var_size_cutoff:
            ax.scatter([], [], marker = '.', c='lightgrey', s=pval_size, label = pval_label, edgecolor = 'grey', linewidth = 0.75)
        else:
            ax.scatter([], [], marker = 'o', c='lightgrey', s=pval_size, label = pval_label, edgecolor = 'k', linewidth = 2)
    plt.legend(scatterpoints=1, frameon=False, labelspacing=1.13, title="pval", markerscale = 1,
               loc = 'lower left', fontsize = 'large', bbox_to_anchor=(1, 0.1), borderaxespad=0)
    ax.tick_params(axis='x', labelrotation=xaxis_rot)
    ax.tick_params(axis='y', labelrotation=yaxis_rot)
#    # place the fold change colorbar
    cbaxes = fig.add_axes([0.93, 0.5, 0.03, 0.3])
    plt.colorbar(sc, shrink= 0.5, aspect = 10, cax = cbaxes)
    #    plt.tight_layout()
    fig.savefig(os.path.join(path, fname), bbox_inches = 'tight')
    return None


def cluster_matrix(matrix, nafil = 1, row_cluster=True, col_cluster=False):
    mask = matrix.isnull()
    cluster_map = sns.clustermap(matrix.fillna(nafil), mask = mask, \
                                 row_cluster=row_cluster, col_cluster=col_cluster, visible=False)
    plt.close()
    if row_cluster == True:
        row_order = cluster_map.dendrogram_row.reordered_ind
        new_index = [matrix.index[i] for i in row_order]
        matrix_reorder = matrix.reindex(new_index)
    if col_cluster == True:
        col_order = cluster_map.dendrogram_col.reordered_ind
        new_col = [matrix.columns[i] for i in col_order]
        if row_cluster != True:
            matrix_reorder = matrix[new_col]
        elif row_cluster == True:
            matrix_reorder = matrix_reorder[new_col]
    if (row_cluster != True) and (col_cluster != True):
        matrix_reorder = matrix.copy(deep=True)
    return matrix_reorder

# function used to assign value to string and return those map and the unique element
def str_to_num(input_str_list, order = 'SortR'):
    if order == 'SortR':
        unique = sorted(set(input_str_list), reverse=True)
        x = [unique.index(str_i) for i, str_i in enumerate(input_str_list)]
    else:
        unique, x = np.unique(input_str_list, return_inverse = True)
        if order == 'AsIsR': 
            x = [max(x)-i for i in x]
    return unique, x


# create dictionary from two list where the first list the key, the second list the value
def dict_twolist(key, value):
    # key and value can both contain duplicate items, but that duplicate item needs to be internally consistent
    out = dict()
    for i, j in zip(key, value):
        if i not in out.keys():
            out[i] = j
    return out


# create a custom colormap where certain cutoff corresponds to different colors
def custom_cmap(st, cutoff_value=[0,0.05,1], cutoff_color=None):
    # cutoff_value and cutoff_color has to be of the same dimension where one cutoff corresponds to one color
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # create the second colorscale called "grey_scale"
    if cutoff_color == None:
        if st == 'AGG':
            colors_1 = [(0.99609375, 0.625, 0.4765625), (0.859375, 0.078125, 0.234375)] ### 0 is 'lightsalmon' ('#FFA07A'), 0.05 is 'crimson'
        elif st == 'TL':
            colors_1 = [(0.390625, 0.58203125, 0.92578125), (0, 0, 0.5)] ### 0 is 'lightskyblue', 0.05 is 'navy'
        elif st == 'Prop':
            colors_1 = [(0.578125, 0, 0.82421875), (0.86328125, 0.625, 0.86328125)] ### 0 is 'darkviolet' ('#9400D3'), 0.05 is 'plum'
        colors_2 = [(0.75, 0.75, 0.75), (1, 1, 1)]
    else:
        colors_1 = cutoff_color[0]
        colors_2 = cutoff_color[1]
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N *(1-cutoff_value[1])/cutoff_value[1])
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)]
    cmap_list += [cmap_2(i) for i in range(cmap_2.N)] 
    cmap_bounds = np.append(np.linspace(0,cutoff_value[1],100), np.linspace(cutoff_value[1],1,100*int((1-cutoff_value[1])/cutoff_value[1])))
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
    return cmap, norm 


################### SECTION 1.1: PLOT BUBBLE PLOT ON CHAPERONINES/Hsp90/Hsp70/Hsp40 ############### 
## old comment: read file of interest, test on chaperonines, the optimal figsize is (6,8) with a size_scale of 180
# updated on 2019.12.01    
# age_comp = 'OvY' # choose among OvY, TvO, TvY
# sample_type = 'TL'
# for set_name in ['Chaperonines', 'Hsp40', 'Hsp70', 'Hsp90', 'NEFs', 'Other_chaperones', 'PDIs', 'PPIases', 'sHsp', 'StressFactors', 'TPRs']:
#     foi_name = '%s_All_%s_%s_flat.csv'%(set_name, age_comp, sample_type)
#     foi_path = 'GeneSets/MS/Chaperones_%s'%sample_type
#     var_value = '%s_logFC'%age_comp; var_size = '%s_pval'%age_comp
#     out_fname = foi_name.replace('flat.csv', 'bubble.pdf')
#     if os.path.isfile(os.path.join(foi_path, foi_name)) == True:
#         df_foi = pd.read_csv(os.path.join(foi_path, foi_name))
#         # create a bubble plot where both x and y went through clustering
#         short_index = 'Human'; short_coln = 'Tissue'
#         df_short = df_foi.pivot_table(index = short_index, columns = short_coln, values = var_value)
#         if (df_short.shape[0]>1) and (df_short.shape[1]>1):
#             df_short_cluster = cluster_matrix(df_short, nafil = 0, row_cluster=True, col_cluster=False)
#             ###### OPTION 1: USE EXISITING TISSUE ORDER
#             #x_unique, x_index = str_to_num(df_short_cluster.columns, order = 'AsIs')
#             ##### OPTION 2: USE SPECIFIC TISSUE ORDER
#             if 'T' in age_comp:
#                 x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Muscle'], [0,1,2,3,4,5]
#             else:
#                 x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Testis', 'Muscle'], [0,1,2,3,4,5,6]
#             x_dict = dict(zip(x_unique, x_index))
#             ##### End of Tissur order specification
#             y_unique, y_index = str_to_num(df_short_cluster.index.tolist(), order = 'AsIs')
#             y_dict = dict(zip(y_unique, y_index))
#             df_long_cluster = df_foi.dropna(subset=([var_value,var_size]))
#             bubble(df_long_cluster, 'Tissue', 'Human', var_value, var_size, \
#                   out_fname.replace('.pdf','_ycluster.pdf'), foi_path, order = 'AsIs',\
#                   var_x_order = x_dict, var_y_order = y_dict, xaxis_rot=90, yaxis_rot = 0)
#            # filter the dataframe to include only the relevant terms
#            bubble_filter(df_long_cluster, 'Tissue', 'Human', var_value, var_size, out_fname.replace('.pdf', '_fil.pdf'),\
#                   foi_path, 'Human', var_x_order = x_dict, order = 'AsIsR', xaxis_rot = 90, yaxis_rot = 0)
############################# END of SECTION 1.1 ############################# 

############## SECTION 1.2: PLOT BUBBLE PLOT ON All Chaperones ############### 
# updated on 2020.04.27    
age1 = 'Young'; age2 = 'Old'
age_comp = '%sv%s'%(age2[0], age1[0])
metric = 'TL' # choose among TL, AGG, and Prop
for set_name in ['Proteasome_All']:#Proteasome_All']:#'AllChaperones']:
    foi_name = '%s_%s_%s_flat.csv'%(set_name, age_comp, metric)
    foi_path = 'GeneSets/MS/%s_%s'%(set_name[0:set_name.index('_')], metric)
    if metric == 'Prop':
        var_value = '%s_prop_logFC'%age_comp; var_size = '%s_prop_pval'%age_comp;
    else:
        var_value = '%s_logFC'%age_comp; var_size = '%s_pval'%age_comp;
    out_fname = foi_name.replace('flat.csv', 'bubble.pdf')
    if os.path.isfile(os.path.join(foi_path, foi_name)) == True:
        df_foi = pd.read_csv(os.path.join(foi_path, foi_name))
        # create a bubble plot where both x and y went through clustering
        short_index = 'Human'; short_coln = 'Tissue'
        df_short = df_foi.pivot_table(index = short_index, columns = short_coln, values = var_value)
        if (df_short.shape[0]>1) and (df_short.shape[1]>1):
            df_short_cluster = cluster_matrix(df_short, nafil = 0, row_cluster=True, col_cluster=False)
            ###### OPTION 1: USE EXISITING TISSUE ORDER
            #x_unique, x_index = str_to_num(df_short_cluster.columns, order = 'AsIs')
            ##### OPTION 2: USE SPECIFIC TISSUE ORDER
            if 'T' in age_comp:
                x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Muscle'], [0,1,2,3,4,5]
            else:
                x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Testis', 'Muscle'], [0,1,2,3,4,5,6]
            x_dict = dict(zip(x_unique, x_index))
            ##### End of Tissur order specification
            y_unique, y_index = str_to_num(df_short_cluster.index.tolist(), order = 'AsIs')
            y_dict = dict(zip(y_unique, y_index))
            df_long_cluster = df_foi.dropna(subset=([var_value,var_size]))
            bubble(df_long_cluster, 'Tissue', 'Human', var_value, var_size, \
                  out_fname.replace('.pdf','_ycluster.pdf'), foi_path, order = 'AsIs',\
                  var_x_order = x_dict, var_y_order = y_dict, xaxis_rot=90, yaxis_rot = 0)
            # filter the dataframe to include only the relevant terms
            if 'All' in set_name:
                bubble_filter(df_long_cluster, 'Tissue', 'Human', var_value, var_size, out_fname.replace('.pdf', '_fil.pdf'),\
                  foi_path, 'Human', var_x_order = x_dict, var_y_order = 'classification', order = 'AsIsR', xaxis_rot = 90, yaxis_rot = 0)
############################# END of SECTION 1.2 ############################# 
                

################# SECTION 2: Plot Bubble Plot On Autophagy Terms ###############
# updated on 2019.12.01    
# age_comp = 'OvY'
# for set_name in ['AllAutophagy']:#, 'AutophagyMitochondrion', 'AutophagyNucleus', 'AutophagyPeroxisome', 'ChaperoneMediatedAutophagy','NegativeRegAutophagy','PositiveRegAutophagy']:
#     foi_name = '%s_%s_AGG_flat.csv'%(set_name, age_comp)
#     foi_path = 'GeneSets/MS/Autophagy_AGG'
#     var_value = '%s_logFC'%age_comp; var_size = '%s_pval'%age_comp
#     out_fname = foi_name.replace('flat.csv', 'bubble.pdf')
#     if os.path.isfile(os.path.join(foi_path, foi_name)) == True:
#         df_foi = pd.read_csv(os.path.join(foi_path, foi_name))
#         # create a bubble plot where both x and y went through clustering
#         short_index = 'Human'; short_coln = 'Tissue'
#         df_short = df_foi.pivot_table(index = short_index, columns = short_coln, values = var_value)
#         if (df_short.shape[0]>1) and (df_short.shape[1]>1):
#             df_short_cluster = cluster_matrix(df_short, nafil = 0, row_cluster=True, col_cluster=False)
#             ###### OPTION 1: USE EXISITING TISSUE ORDER
#             #x_unique, x_index = str_to_num(df_short_cluster.columns, order = 'AsIs')
#             ##### OPTION 2: USE SPECIFIC TISSUE ORDER
#             if 'T' in age_comp:
#                 x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Muscle'], [0,1,2,3,4,5]
#             else:
#                 x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Testis', 'Muscle'], [0,1,2,3,4,5,6]
#             x_dict = dict(zip(x_unique, x_index))
#             ##### End of Tissur order specification
#             y_unique, y_index = str_to_num(df_short_cluster.index.tolist(), order = 'AsIs')
#             y_dict = dict(zip(y_unique, y_index))
#             df_long_cluster = df_foi.dropna(subset=([var_value,var_size]))
#             bubble(df_long_cluster, 'Tissue', 'Human', var_value, var_size, \
#                   out_fname.replace('.pdf','_ycluster.pdf'), foi_path, order = 'AsIs',\
#                   var_x_order = x_dict, var_y_order = y_dict, xaxis_rot=90, yaxis_rot = 0)
#             # filter the dataframe to include only the relevant terms
#             if 'All' in set_name:
#                 bubble_filter(df_long_cluster, 'Tissue', 'Human', var_value, var_size, out_fname.replace('.pdf', '_fil.pdf'),\
#                   foi_path, 'Human', var_x_order = x_dict, var_y_order = 'classification', order = 'AsIsR', xaxis_rot = 90, yaxis_rot = 0)   
################################# END of SECTION 2 #############################  


################ SECTION 3: Do clustering on the results, then plot bubble plot on autoghagy terms ############### 
## update on 2019.12.01
# age_comp = 'OvY'
# set_name = 'NegativeRegAutophagy' #PositiveRegAutophagy
# foi_name = '%s_%s_AGG_flat.csv'%(set_name, age_comp)
# foi_path = 'GeneSets/MS/%s_AGG'%set_name.replace('NegativeReg','').replace('PositiveReg','')
# var_value = '%s_logFC'%age_comp; var_size = '%s_pval'%age_comp
# out_fname = foi_name.replace('flat.csv', 'bubble.pdf').replace('_All', '')
# df_foi = pd.read_csv(os.path.join(foi_path, foi_name))
# ## create a normal bubble plot
# #bubble(df_foi.dropna(subset=([var_value,var_size])), 'Tissue', 'Human', var_value, var_size, \
# #       out_fname.replace('.pdf','_original.pdf'), foi_path, order = 'AsIs', xaxis_rot=90, yaxis_rot = 0)
# # create a bubble plot where both x and y went through clustering
# short_index = 'Human'; short_coln = 'Tissue'
# df_short = df_foi.pivot_table(index = short_index, columns = short_coln, values = var_value)
# df_short_cluster = cluster_matrix(df_short, nafil = 0, row_cluster=True, col_cluster=False)
# ###### OPTION 1: USE EXISITING TISSUE ORDER
# #x_unique, x_index = str_to_num(df_short_cluster.columns, order = 'AsIs')
# ##### OPTION 2: USE SPECIFIC TISSUE ORDER
# if 'T' in age_comp:
#     x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Muscle'], [0,1,2,3,4,5]
# else:
#     x_unique, x_index = ['Heart', 'Brain', 'Liver', 'Gut', 'Skin', 'Testis', 'Muscle'], [0,1,2,3,4,5,6]
# x_dict = dict(zip(x_unique, x_index))
# ##### End of Tissur order specification
# y_unique, y_index = str_to_num(df_short_cluster.index.tolist(), order = 'AsIs')
# y_dict = dict(zip(y_unique, y_index))
# df_long_cluster = df_foi.dropna(subset=([var_value,var_size]))
# bubble(df_long_cluster, 'Tissue', 'Human', var_value, var_size, \
#       out_fname.replace('.pdf','_ycluster.pdf'), foi_path, order = 'AsIs',\
#       var_x_order = x_dict, var_y_order = y_dict, xaxis_rot=90, yaxis_rot = 0)
# bubble_filter(df_long_cluster, 'Tissue', 'Human', var_value, var_size, out_fname.replace('.pdf', '_ycluster_fil_ver.pdf'),\
#       foi_path, 'Human', var_x_order = x_dict, xaxis_rot = 0, yaxis_rot = 0)
############################# END of SECTION 5 #############################