#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 23:48:31 2019

this is used to generate heatmap for the correlomics/enrichment differences among different parameters

@author: yiwenchen
"""

from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import itertools
import sys


mpl.rcParams['pdf.fonttype'] = 42
#sns.set(style="white")


# this version is when the darker the color, the more significant the p-value 
def two_scale_heatmap(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = True, map_square = True, \
                      figsize = (10,10), yaxis_rot = 0):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    f, ax = plt.subplots(figsize=figsize)
    if st == 'AGG':
        colors_1 = [(0.99609375, 0.625, 0.4765625), (0.859375, 0.078125, 0.234375)] ### 0 is 'lightsalmon' ('#FFA07A'), 0.05 is 'crimson'
    elif st == 'TL':
        colors_1 = [(0.390625, 0.58203125, 0.92578125), (0, 0, 0.5)] ### 0 is 'lightskyblue', 0.05 is 'navy'
    else:
        colors_1 = [(0.578125, 0, 0.82421875), (0.86328125, 0.625, 0.86328125)] ### 0 is 'darkviolet' ('#9400D3'), 0.05 is 'plum'
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    if colors_2 == 'None': colors_2 = [(0.75, 0.75, 0.75), (1, 1, 1)]
    # create the second colorscale called "grey_scale"
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N *(1-cutoff)/cutoff)
    #cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N * 9)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)]
    cmap_list += [cmap_2(i) for i in range(cmap_2.N)] 
    cmap_bounds = np.append(np.linspace(0,cutoff,100), np.linspace(cutoff,1,100*int((1-cutoff)/cutoff)))
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    sns.set(style="white")
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
                square=map_square, ax = ax, linewidths=.5, cbar_kws={"shrink": .5}, yticklabels = ylabels) 
    # without linewidth=.5, the figure looks much better
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    plt.close()


# this version is when there are two kinds of shade for the heatmap, each representing a different change
# in the direction of the variable (i.e. increase with green, decrease with pink determined by t stats), 
# then the actual shade of the individual grid is indicative of a variable value (i.e the p-value) 
def two_direction_heatmap(matrix, out_path, out_fname, mask='None', st = 'None', \
                          colors = None, cbar=True, positions = None, ylabels = True, map_square = True, \
                      figsize = (10,10), yaxis_rot = 0, vlim = [-1,1], title = None, annot=False, ax = None):
    if (colors == None) & (positions == None):
        if 'pval' in  st:
            cmap = make_cmap([(1,1,1), (0.75,0.75,0.75), (0.390625, 0.58203125, 0.92578125),\
                                 (0, 0, 0.5),(0.859375, 0.078125, 0.234375),\
                                 (0.99609375, 0.625, 0.4765625), (0.75,0.75,0.75), (1,1,1)], \
                        position=[0,0.4749999,0.475,0.49999,0.5,0.5249999,0.525,1], bit=False)
        elif 'zscore' in st:
#            cmap = make_cmap([(0, 0, 1),(1, 1, 1),(1, 0, 0)], position=[0,0.5,1], bit=False)  # this is blue-white-red
            cmap = make_cmap([(0, 0, 1),(1, 1, 1),(1, 1, 0)], position=[0,0.5,1], bit=False)  # this is blue-white-yellow
    else:
        cmap = colors
    mask = np.zeros_like(matrix, dtype=np.bool) if mask == 'None' else matrix.isnull() # if True is no mask
    if ax == None: 
        f, ax = plt.subplots(figsize=figsize)
        out = True
    else:
        out = False
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 3)
    # detailed style guide: https://seaborn.pydata.org/tutorial/aesthetics.html
    sns.set_style("darkgrid", {'axes.facecolor':'lightgrey'})
    matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin = vlim[0], vmax=vlim[1], annot=annot, fmt='.2g', annot_kws={"size": 10},\
                square=map_square, ax = ax, linewidths=2, cbar = cbar, cbar_kws={"shrink": .5}, yticklabels = ylabels) 
    if title != None: ax.set(title = title)
    # without linewidth=.5, the figure looks much better
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    if out == True: matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
    plt.close()


# this version is when the lighter the color, the more significant the p-value 
def two_scale_heatmap_light(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0):
    mask = np.zeros_like(matrix, dtype=np.bool) if mask == 'None' else matrix.isnull() # if True is no mask
    f, ax = plt.subplots(figsize=figsize)
    cmap_1 = plt.cm.Purples
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html
    colors_2 = [(0.66, 0.66, 0.66), (0.3, 0.3, 0.3)]
    # create the second colorscale called "grey_scale"
    cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N * 9)
    #cmap_2 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_2, N = cmap_1.N * 9)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)]
    cmap_list += [cmap_2(i) for i in range(cmap_2.N)] 
    cmap_bounds = np.append(np.linspace(0,0.1,100), np.linspace(0.1,1,900))
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
    ####
    # Draw the heatmap with the mask and correct aspect ratio
    sns.set(font_scale = 4)
    sns.heatmap(matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
                square=True, linewidths=.5, ax = ax, cbar_kws={"shrink": .5})
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    f.savefig(os.path.join(out_path, out_fname.replace('.csv', '_heatmap.pdf')))
    plt.close(f)


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


def heatmap_hybrid(matrix, out_path, out_fname, cutoff = 0.05, mask='None', st = 'None', \
                      colors_1 = 'None', colors_2 = 'None', ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0, row_cluster=True, col_cluster=False):
    matrix_reorder = cluster_matrix(matrix, row_cluster=row_cluster, col_cluster=col_cluster)
    sns.set(font_scale = 2)
    two_scale_heatmap(matrix_reorder, out_path, out_fname, cutoff = cutoff, mask=mask, st = st, \
                      colors_1 = colors_1, colors_2 = colors_2, ylabels = ylabels, \
                      map_square = map_square, figsize = figsize, yaxis_rot = yaxis_rot) 

def heatmap2d_hybrid(matrix, out_path, out_fname, mask='None', st = 'None', \
                      colors = None, positions = None, ylabels = False, map_square = False, \
                      figsize = (10,10), yaxis_rot = 0, row_cluster=True, col_cluster=False):
    matrix_reorder = cluster_matrix(matrix, row_cluster=row_cluster, col_cluster=col_cluster)
    sns.set(font_scale = 2)
    two_direction_heatmap(matrix_reorder, out_path, out_fname, mask=mask, st = st, colors = colors,\
                      positions = positions, ylabels = ylabels, map_square = map_square, \
                      figsize = figsize, yaxis_rot = yaxis_rot)


def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))
    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap


# Generate a custom diverging colormap
#cmap = sns.diverging_palette(220, 10, as_cmap=True)
#cmap = sns.diverging_palette(360, 10, as_cmap = True, l=0, s=10)


####### OPTION 1: Generate a custom heatmap where it is linear up till a cutoff, then grey out for the rest ######
#cmap = plt.cm.BuPu
#cmap_list = [cmap(i) for i in range(cmap.N)]
#cmap_list += [(0.5, 0.5, 0.5, 1.0)]*cmap.N*9
#cmap_bounds = np.append(np.linspace(0,0.1,101), np.linspace(0.1,1,901))
## create the new colormap
#cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap.N)
#norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
################################################# END OF OPTION 1 ################################################


####### OPTION 2: Generate a custom heatmap where it is linear up till a cutoff, then grey scale for the rest ######
#cmap_1 = plt.cm.Blues
##colors_2 = [(0.5, 0.5, 0.5), (0, 0, 0)]
#cmap_2 = plt.cm.Greys
#cmap_list = [cmap_1(i) for i in range(cmap_1.N)]
#for i in range(cmap_2.N):
#    cmap_list += [cmap_2(i)]*9
#cmap_bounds = np.append(np.linspace(0,0.1,100), np.linspace(0.1,1,900))
## create the new colormap
#cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
#norm = mpl.colors.BoundaryNorm(cmap_bounds, cmap.N)
################################################ END OF OPTION 2 ################################################


########## Section 3: Generate heatmap for biophysical properties that went through variable stats test ##########
## direction of changes as determined by the t-test stats is used to generate two colorscale 
### part 2: plot the permutation results update on 2020.04.26 (focus on age-comparison)
weights = 'Young_count'  #  None or 'count_Young_TL' or 'Young_count'
weights_type = '' if weights is None else '_%s_weights'%weights
query_metrics = ['AGG','Prop'] if weights is None else ['AGG']
p_cutoff =  .05
age_comp = 'OvY' # OvY or TvO
test_type = 'ks-binary'#'ks'
map_value = 'hit.mean.zscore';  map_mask = 'permutation.pval'
for metric in query_metrics:
    # read the hit dataframe
    hit_fpath = 'Correlomics/%sHits_%s/Permutation'%(metric, age_comp) if weights is None else 'Correlomics/%sHits_%s/PermutationWeights'%(metric, age_comp) 
    for hit_dir in ['Pos', 'PosZ0.75']:#['Neg', 'NegZ0.75']:#, 'Pos', 'PosZ0.75']:
        fname_hit = '%sSig%s_%s_%s%s'%(metric, hit_dir, age_comp, test_type, weights_type)
        df_hit = pd.read_csv(os.path.join(hit_fpath, '%s.csv'%fname_hit))
        df_wide = pd.pivot_table(df_hit, index='Tissue', columns = 'Metric', values = [map_value, map_mask]) 
        df_mask = df_wide[map_mask]
        df_result = df_wide[map_value]
#        ##### Option 1 where the color is based on p-value and the direction of  change is reflected in the color
#        # result_adj is made such that a negative value means it is decreased with age, still plot based on p-value
#        df_result_adj = pd.DataFrame(data = np.where(df_direction.values<0, df_result.values*-1, df_result.values),\
#                                     index = df_result.index, columns = df_result.columns)
        ##### Option 2 where the color is based on the zscore
        df_result_adj = pd.DataFrame(data = np.where(df_mask.values<=p_cutoff, df_result.values, np.nan),\
                                    index = df_result.index, columns = df_result.columns)
        print ('max is %s, min is %s'%(df_result_adj.max().max(), df_result_adj.min().min()))
        # make the relevant heatmap
#        col_order = ['Frac_pos', 'FracNeg','FracCharged','FracExpanding','delta',\
#                       'FracDisoPromoting', 'Disordered Frac', 'AA Length',\
#                       'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu','Frac_QN','VITmaxrun', 'PRDscore','MeanNetCharge',\
#                       'mean Ka/Ks','mean Ka','mean ks',\
#                       'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
#                       'log2_Young_TL','log2_Old_TL','log2Prop_Young','OvY_logFC', 'Delta_Propensity_Old-Young']
        col_order = ['log2_Young_TL', 'AA Length', 'FImaxrun','Longest Disorder Stretch',\
                    'FracExpanding','FracDisoPromoting', 'Disordered Frac','Frac_FInumaa',\
                    'PAPAprop','NLLR','MWScore','Frac_QN',\
                    'pI', 'FracNeg','Frac_pos','FracCharged',\
                    'delta', 'kappa','Omega','MeanNetCharge','uversky_hydropathy','Frac_neu',\
                    'Frac_aromatic','Frac_aliphatic',\
                    'mean Ka/Ks','mean Ka','mean ks']# for the supp figure,'OvY_logFC', 'Delta_Propensity_Old-Young']
#        col_order = ['NLLR','Disordered Frac']
        if 'T' in age_comp:
#            row_order = ['Skin','Brain','Gut','Liver','Heart','Muscle']
            row_order = ['Skin','Muscle','Liver','Heart','Gut','Brain']
        else:
#            row_order = ['Skin','Brain','Gut','Liver','Testis','Heart','Muscle']
#            row_order = ['Heart','Brain','Liver','Gut','Skin','Testis','Muscle']
            row_order = ['Testis','Skin','Muscle','Liver','Heart','Gut','Brain']
#        sns.set_style("darkgrid", {'axes.facecolor':'lightgrey'})
        two_direction_heatmap(df_result_adj.loc[row_order, col_order], hit_fpath, '%s_p%s_heatmap2d_custom.pdf'%(fname_hit, p_cutoff), st = map_value, \
                  figsize = (24, 15), yaxis_rot = 0, map_square = True, vlim = [-3,3], ylabels=True)
##### part 3: plot the permutation results update on 2020.05.01 (focus on abundant AGG and/or Prop)
weights = None  #  None or 'count_Young_TL'
weights_type = '' if weights is None else '_%s_weights'%weights
query_metrics = ['AGG','Prop'] if weights is None else ['AGG']
p_cutoff =  .05; correction = False
z_cutoff = 2
age = 'Young'
test_type = 'ks-binary'#'ks'
map_value = 'hit.mean.zscore';  map_mask = 'permutation.pval'
for metric in query_metrics:
    # read the hit dataframe
    hit_fpath = 'Correlomics/AgeSpecific/Permutation' if weights is None else 'Correlomics/AgeSpecific/PermutationWeights' 
    for hit_dir in ['Pos', 'PosZ0.75']:#['Neg', 'NegZ0.75']:#, 'Pos', 'PosZ0.75']:
        fname_hit = '%s_%sZ%s_%s%s'%(age, metric, z_cutoff, test_type, weights_type)
        df_hit = pd.read_csv(os.path.join(hit_fpath, '%s.csv'%fname_hit))
        df_wide = pd.pivot_table(df_hit, index='Tissue', columns = 'Metric', values = [map_value, map_mask]) 
        df_mask = df_wide[map_mask]
        df_result = df_wide[map_value]
#        ##### Option 1 where the color is based on p-value and the direction of  change is reflected in the color
#        # result_adj is made such that a negative value means it is decreased with age, still plot based on p-value
#        df_result_adj = pd.DataFrame(data = np.where(df_direction.values<0, df_result.values*-1, df_result.values),\
#                                     index = df_result.index, columns = df_result.columns)
        ##### Option 2 where the color is based on the zscore
        col_order = ['log2_Young_TL', 'AA Length', 'FImaxrun','Longest Disorder Stretch',\
                     'FracExpanding','FracDisoPromoting', 'Disordered Frac','Frac_FInumaa',\
                     'PAPAprop','NLLR','MWScore','Frac_QN',\
                     'pI', 'FracNeg','Frac_pos','FracCharged',\
                     'delta', 'kappa','Omega','MeanNetCharge','uversky_hydropathy','Frac_neu',\
                     'Frac_aromatic','Frac_aliphatic',\
                     'mean Ka/Ks','mean Ka','mean ks']# for the supp figure,'OvY_logFC', 'Delta_Propensity_Old-Young']
        if 'T' in age:
            row_order = ['Skin','Muscle','Liver','Heart','Gut','Brain']
        else:
            row_order = ['Testis','Skin','Muscle','Liver','Heart','Gut','Brain']
        if correction == True: p_cutoff = p_cutoff/len(row_order)/len(col_order)
        df_result_adj = pd.DataFrame(data = np.where(df_mask.values<=p_cutoff, df_result.values, np.nan),\
                                     index = df_result.index, columns = df_result.columns)
        print ('max is %s, min is %s'%(df_result_adj.max().max(), df_result_adj.min().min()))
        # make the relevant heatmap
#        col_order = ['Frac_pos', 'FracNeg','FracCharged','FracExpanding','delta',\
#                       'FracDisoPromoting', 'Disordered Frac', 'AA Length',\
#                       'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu','Frac_QN','VITmaxrun', 'PRDscore','MeanNetCharge',\
#                       'mean Ka/Ks','mean Ka','mean ks',\
#                       'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
#                       'log2_Young_TL','log2_Old_TL','log2Prop_Young','OvY_logFC', 'Delta_Propensity_Old-Young']
#        sns.set_style("darkgrid", {'axes.facecolor':'lightgrey'})
        two_direction_heatmap(df_result_adj.loc[row_order, col_order], hit_fpath, '%s_p%s_heatmap2d_custom.pdf'%(fname_hit, p_cutoff), st = map_value, \
                  figsize = (24, 15), yaxis_rot = 0, map_square = True, vlim = [-8,8], ylabels=True)
######### End of Section 5 ##########
########## End of Section 5 ##########