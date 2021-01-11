import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import numpy as np
import math
import matplotlib as mpl
import scipy
from scipy import stats
import random
import pickle
import time

mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 2)
sns.set(style="white")



def make_kde(*args, **kwargs):    
    sns.kdeplot(*args, cmap=next(make_kde.cmap_cycle), **kwargs)    


def spearmanfunc(x, y, **kwargs):
    rho, pval = stats.spearmanr(x, y)
    ax = plt.gca()
    ax.annotate('r = %.2f, p-value = %.2f'%(rho, pval), xy=(.1, .9), xycoords = ax.transAxes)   


def hist_samebin(metric, df_data, ax, hue_coln, hue_order = None, colors_list = None,\
                 alpha = 1, histtype = 'step', bins= 10, lw = 1, normsum = False, weight_coln = None, ks_test = False, legend=True):
    df_plot = df_data[metric] if weight_coln == None else df_data[metric]*df_data[weight_coln]
    metric_min = df_plot.min()#df_plot[metric].min()
    metric_max = df_plot.max()
    for hue_i, color in zip(hue_order, colors_list):
        hist_x = df_plot[df_data[hue_coln]==hue_i].dropna()
        if normsum == False:
            # # use histplot for python 3., uncomment the next line if using python 2. because histplot was introduced much later
            # sns.histplot(hist_x.values, kde = False, ax = ax, bins = bins, stat ='count', fill = False, \
            #           element = histtype, color = color, alpha = alpha, line_kws={'lw':lw},\
            #                         label = hue_i, binrange = (metric_min, metric_max))
            # use distplot for python 2, this is depreciated in python 3 for later versions of seaborn    
            sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins = bins,\
                      hist_kws = {'histtype':histtype,'alpha':alpha,'lw':lw,'color':color,\
                                  'range':(metric_min, metric_max),'label':hue_i})#,'density': True})
        # normsum is to normalize the sum of the bin counts to total of 1 (i.e. histogram on fraction of counts)
        elif normsum == True:
            # # use histplot for python 3., uncomment the next line if using python 2. because histplot was introduced much later
            # sns.histplot(hist_x.values, kde = False, ax = ax, bins = bins, stat = 'probability', fill = False, \
            #                         element = histtype, color = color, alpha = alpha, line_kws={'lw':lw},\
            #                         label = hue_i, binrange = (metric_min, metric_max))   
            # use distplot for python 2, this is depreciated in python 3 for later versions of seaborn
            sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins=bins,\
                          hist_kws = {'histtype':histtype, 'alpha':alpha, 'lw':lw, 'color':color,\
                                      'label':hue_i, 'weights': np.ones_like(hist_x)*1./len(hist_x), 'range':(metric_min, metric_max)})
    if ks_test == True:
            if df_data[hue_coln].nunique() == 2:
                hue_col_vs = df_data[hue_coln].unique().tolist()
                ks_stat, ks_pvalue = stats.ks_2samp(df_data[metric][df_data[hue_coln]==hue_col_vs[0]],df_data[metric][df_data[hue_coln]==hue_col_vs[1]])
#                fig.text(0.5, 0.9, 'ks: %.2f, p-value: %.2f'%(ks_stat, ks_pvalue), horizontalalignment='center', size='medium', color='black', transform=ax.transAxes)
                ax.set_title('%s\n ks = %.2f, p-value = %.2f'%(metric, ks_stat, ks_pvalue))
    if legend == True: ax.legend()
    return ax


def hist_cum(metric, df_data, ax, hue_coln, hue_order = None, colors_list = None,\
                 alpha = 1, bins= 10, lw = 1, normsum = False, weight_coln = None, ks_test = False):
    df_plot = df_data[metric] if weight_coln == None else df_data[metric]*df_data[weight_coln]
    metric_min = df_plot.min()#df_plot[metric].min()
    metric_max = df_plot.max()
    for hue_i, color in zip(hue_order, colors_list):
        hist_x = df_plot[df_data[hue_coln]==hue_i].dropna()
        # use histplot for python 3., uncomment the next line if using python 2. because histplot was introduced much later
        sns.histplot(hist_x.values, kde = False, ax = ax, bins = bins, cumulative=True, \
                      element = 'step', stat = 'density', alpha = alpha, fill = False, \
                      line_kws={'lw':lw},color=color, binrange = (metric_min, metric_max), label = hue_i)
        # # use distplot for python 2, this is depreciated in python 3 for later versions of seaborn
        # sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins = bins, \
        #              hist_kws = {'histtype':'step', 'cumulative': True,'alpha':alpha,'lw':lw,'color':color,\
        #                          'range':(metric_min, metric_max),'label':hue_i, 'density': True})
    if ks_test == True:
            if df_data[hue_coln].nunique() == 2:
                hue_col_vs = df_data[hue_coln].unique().tolist()
                ks_stat, ks_pvalue = stats.ks_2samp(df_data[metric][df_data[hue_coln]==hue_col_vs[0]],df_data[metric][df_data[hue_coln]==hue_col_vs[1]])
#                fig.text(0.5, 0.9, 'ks: %.2f, p-value: %.2f'%(ks_stat, ks_pvalue), horizontalalignment='center', size='medium', color='black', transform=ax.transAxes)
                ax.set_title('%s\n ks = %.2f, p-value = %.2f'%(metric, ks_stat, ks_pvalue))
    return ax


# make pair-wise scatter plot 
def grid_correlomics(fig_name, df_data, matrix_coln, info_coln, hue_coln = None, fig_path ='', hue_order = False, colors_list = False, alpha = 1):
    colns = matrix_coln + info_coln 
    print (df_data.shape)
    df_plot = df_data[colns].dropna()
    print (df_plot.shape)
    if colors_list == False: 
        colors = sns.color_palette('husl', 10)
        colors_list = colors.as_hex()
    sns.set(style = "ticks")
    if hue_coln == None:
        g = sns.PairGrid(df_plot, palette = sns.color_palette(colors_list))
    else:
        if hue_order == False:
            g = sns.PairGrid(df_plot, hue = hue_coln, palette = sns.color_palette(colors_list), diag_sharey= False)
        else:
            g = sns.PairGrid(df_plot, hue = hue_coln, hue_order = hue_order, palette = sns.color_palette(colors_list),diag_sharey= False)
    # plot histogram in the diagonal space
    g = g.map_diag(plt.hist, bins=10, histtype = 'step', linewidth = 3)#, density = True) 
#    for i, coln in enumerate(matrix_coln):
#        hist_samebin(coln, df_data, g.axes[i,i], hue_coln, hue_order = hue_order, colors_list = colors_list, normsum = True)
#        g.axes[i,i].set_ylim((0,1))
    # plot scatter plot on the upper corner
    g = g.map_upper(plt.scatter, alpha = alpha)
    make_kde.cmap_cycle = itertools.cycle([sns.light_palette(i, as_cmap=True) for i in colors_list])
    # plot kde in the lower corner (display spearman coefficient if appropriate)
    g = g.map_lower(make_kde)
    if hue_coln == None:
        g.map_lower(spearmanfunc)
    elif df_data[hue_coln].nunique() == 1:
        g.map_lower(spearmanfunc)
    g = g.add_legend() 
    g.savefig(os.path.join(fig_path, '%s.pdf'%fig_name))
    plt.close()    


# make histogram from dataframe with hue
def hist_hue(df_data, hist_colns, hist_fname, hue_coln, hue_order = None, colors_list = None, hit_count=None, alpha = 1, 
             sharey = False, histtype = 'step', bins= 10, lw = 3, normsum = False, ks_test = False, directory = ''):
    if hue_order == None:
        hue_order = sorted(set(df_data[hue_coln].values))
        if 'null' in hue_order: hue_order.append(hue_order.pop(hue_order.index('null')))
    if (hue_order != None) & (colors_list != None):
        if len(hue_order) != len(colors_list):
            print ('Default colors used instead because the specified color do not match with hue_coln')
            colors_list = sns.color_palette('husl', len(hue_order))
    if colors_list == None:
        colors_list = sns.color_palette('husl', len(hue_order))
    hist_coln_num = 4;     
    hist_row_num = int(math.ceil(len(hist_colns)/1./hist_coln_num))
    fig_hist, ax_hist = plt.subplots(nrows = hist_row_num, ncols = hist_coln_num, figsize = (16,4*hist_row_num), sharey = sharey)
    ks_test_out = [[],[], []]
    for ax, metric in zip(ax_hist.flat, hist_colns):
        df_plot = df_data[[metric,hue_coln]].dropna(subset=[metric]).apply(pd.to_numeric,errors='coerce')
        metric_min = df_plot[metric].min()
        metric_max = df_plot[metric].max()
        for hue_i, color in zip(hue_order, colors_list):
            hist_x = df_plot[metric][df_data[hue_coln]==hue_i]
            if hit_count == hue_i: n_hit = len(hist_x)
            if normsum == False:
                # # use histplot for python 3., uncomment the next line if using python 2. because histplot was introduced much later
                # sns.histplot(hist_x.values, kde = False, ax = ax, bins = bins,stat ='count', fill=False, \
                #                    element = histtype, color = color, alpha = alpha, line_kws={'lw':lw},\
                #                    label = hue_i, binrange = (metric_min, metric_max))
                # use distplot for python 2, this is depreciated in python 3 for later versions of seaborn
                sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins = bins,\
                         hist_kws = {'histtype':histtype,'alpha':alpha,'lw':lw,'color':color,\
                                     'range':(metric_min, metric_max),'label':hue_i})
            # normsum is to normalize the sum of the bin counts to total of 1 (i.e. histogram on fraction of counts)
            elif normsum == True:
                # # use histplot for python 3., uncomment the next line if using python 2. because histplot was introduced much later
                # sns.histplot(hist_x.values, kde = False, ax = ax, bins = bins, stat = 'probability', fill=False, \
                #                    element = histtype, color = color, alpha = alpha, line_kws={'lw':lw},\
                #                    label = hue_i, binrange = (metric_min, metric_max))
                # use distplot for python 2, this is depreciated in python 3 for later versions of seaborn
                sns.distplot(hist_x.values, hist = True, kde = False, ax = ax, bins=bins,\
                             hist_kws = {'histtype':histtype,'alpha':alpha,'lw':lw,'color':color,\
                                         'label':hue_i, 'weights': np.ones_like(hist_x)*1./len(hist_x), 'range':(metric_min, metric_max)})
        if ks_test == True:
            if df_data[hue_coln].nunique() == 2:
                hue_col_vs = df_data[hue_coln].unique().tolist()
                ks_stat, ks_pvalue = stats.ks_2samp(df_data[metric][df_data[hue_coln]==hue_col_vs[0]],df_data[metric][df_data[hue_coln]==hue_col_vs[1]])
#                fig.text(0.5, 0.9, 'ks: %.2f, p-value: %.2f'%(ks_stat, ks_pvalue), horizontalalignment='center', size='medium', color='black', transform=ax.transAxes)
                if hit_count!=None:
                    ax.set_title('%s\n ks = %.2f, p-value = %.2f\nn_POI = %s'%(metric, ks_stat, ks_pvalue, n_hit))
                else:
                    ax.set_title('%s\n ks = %.2f, p-value = %.2f'%(metric, ks_stat, ks_pvalue))
                ks_test_out[0].append(metric)
                ks_test_out[1].append(ks_stat)
                ks_test_out[2].append(ks_pvalue)
        else:
            ax.set_title(metric)
    handles,labels = ax.get_legend_handles_labels()
    fig_hist.legend(handles, labels, loc = 'center right')
    plt.tight_layout()
    fig_hist.savefig(os.path.join(directory,'%s.pdf'%hist_fname))
    plt.close(fig_hist)
    return ks_test_out


# make pair-wise scatter plot 
def grid_correlomics_diakde(fig_name, df_data, matrix_coln, info_coln, hue_coln = None, fig_path ='', hue_order = False, colors_list = False, alpha = 1):
    colns = matrix_coln + info_coln 
    df_plot = df_data[colns].dropna()
    if colors_list == False: 
        colors = sns.color_palette('husl', 10)
        colors_list = colors.as_hex()
    sns.set(style = "ticks")
    if hue_coln == None:
        g = sns.PairGrid(df_plot, palette = sns.color_palette(colors_list))
    else:
        if hue_order == False:
            g = sns.PairGrid(df_plot, hue = hue_coln, palette = sns.color_palette(colors_list), diag_sharey= False)
        else:
            g = sns.PairGrid(df_plot, hue = hue_coln, hue_order = hue_order, palette = sns.color_palette(colors_list), diag_sharey= False)
    # plot kde estimated histogram in the diagonal space
    g = g.map_diag(sns.kdeplot, lw=3, legend=True, alpha = alpha)
    # plot scatter plot on the upper corner
    g = g.map_upper(plt.scatter, alpha = alpha)
    make_kde.cmap_cycle = itertools.cycle([sns.light_palette(i, as_cmap=True) for i in colors_list])
    # plot kde in the lower corner (display spearman coefficient if appropriate)
    g = g.map_lower(make_kde)
    if hue_coln == None:
        g.map_lower(spearmanfunc)
    elif df_data[hue_coln].nunique() == 1:
        g.map_lower(spearmanfunc)
    g = g.add_legend() 
    g.savefig(os.path.join(fig_path, fig_name+'.pdf')) 
    plt.close(g)


def hits_all_metric(metric, comp, hit_path, f_path='Properties', z_cutoff = None):
    fname_hits = '%sSig_%s_zcutoff.csv'%(metric, comp) if z_cutoff != None else '%sSig_%s.csv'%(metric, comp)
    fname_all = '%sAll_%s.csv'%(metric, comp)        
    # import the hits list (coln 1 XP_Protein, coln 2 Final_Symbol, coln 3 Tissue, coln 4 Category)
    df_hits = pd.read_csv(os.path.join(hit_path, fname_hits))
    # import the protein list (coln 1 XP_Protein, coln 2 Final_Symbol, coln 3 Tissue)
    df_all = pd.read_csv(os.path.join(hit_path, fname_all))
    # import all the biophysical properties
    df_metrics = pd.read_csv(os.path.join(f_path, 'Killifish_DISOLLRCider.csv'))
    # merge all the biophysical properties with "hit" and "all" table
    df_all_metrics = pd.merge(df_all, df_metrics, how = 'left', left_on = 'N. furzeri Protein Id', right_on = 'N. furzeri Protein Id')
    df_hits_metrics = pd.merge(df_hits, df_metrics, how = 'left', left_on = 'N. furzeri Protein Id', right_on = 'N. furzeri Protein Id')
    # concatenate hits and all matrix
    df_matrix = pd.concat([df_hits, df_all], axis = 0, ignore_index = False, sort = False)
    return df_matrix, df_all_metrics, df_hits_metrics    


def regeneration(tissue):
    if tissue in ['Liver', 'Gut', 'Skin', 'Testis','Muscle']:
        out = 'mitotic'
    elif tissue in ['Brain', 'Heart']:#, 'muscle']:
        out = 'post-mitotic'
    return out


###### this generates a regular heatmap that go from 0 to 1 in one directional color scale ######
def regular_heatmap(matrix, out_path, out_fname, mask='None', st = 'None', annot=False,\
                      colors_1 = 'None', colors_2 = 'None', xlabels = False, ylabels = False, map_square = True, \
                      figsize = (10,10), xaxis_rot = 0, yaxis_rot = 0, vmin=0, vmax = 1, ax = 'None'):
    if mask == 'None': mask = np.zeros_like(matrix, dtype=np.bool) # this is no mask
    if ax == 'None':
        f, ax = plt.subplots(figsize=figsize)
    # (0.5, 0.5, 0.5) is grey in RGBA color format, (0, 0, 0) is black, so start with grey and end with black
    # avoid using the alpha transparency metric because it will screw up the colorbar
    # to find out the rgb value of a color, use webcolors.name_to_rgb('darkgray'), to convert to rgba, divive by 256
    # for list of common color and their string name, visit https://matplotlib.org/examples/color/named_colors.html and 
    # https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    # this version is when the darker the color, the more significant the p-value 
    if st == 'AGG':
        colors_1 = [(1, 1, 1), (0.859375, 0.078125, 0.234375)] ### (0.99609375, 0.625, 0.4765625) is 'lightsalmon', (0.859375, 0.078125, 0.234375) is 'crimson'
    elif st == 'TL':
        colors_1 = [(1, 1, 1), (0, 0, 0.5)] ### (0.390625, 0.58203125, 0.92578125) is 'lightskyblue', (0, 0, 0.5) is 'navy'
    else:
        colors_1 = [(1, 1, 1), (0.578125, 0, 0.82421875)] ### (0.86328125, 0.625, 0.86328125) is 'darkviolet', (0.578125, 0, 0.82421875) is 'plum'
    # create the second colorscale called "grey_scale"
    cmap_1 = mpl.colors.LinearSegmentedColormap.from_list('grey_scale', colors_1, N = 256)
    cmap_list = [cmap_1(i) for i in range(cmap_1.N)] 
    # create the new colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('pval_cmap', cmap_list, cmap_1.N)
    # Draw the heatmap with the mask and correct aspect ratio
#    sns.set(font_scale = 3)
    if annot == False:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, \
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, xticklabels = xlabels, yticklabels = ylabels) 
    else:
        matrix_map = sns.heatmap(matrix, cmap = cmap, mask=mask, vmin=vmin, vmax = vmax, annot=True, fmt='.2g',\
                square=map_square, ax = ax, cbar_kws={"shrink": .5}, xticklabels = xlabels, yticklabels = ylabels)
    #matrix_map = sns.heatmap(df_matrix, cmap = cmap, mask=mask, norm = norm, vmin = 0, vmax=1,\
    #            square=False, linewidths=.5, ax = ax, cbar_kws={"shrink": .5, "orientation": 'horizontal'}, yticklabels = False)
    matrix_map.tick_params(labelsize=24, labelrotation=xaxis_rot)
    matrix_map.tick_params(axis='y', labelrotation=yaxis_rot)
    #sns.clustermap(matrix.dropna(axis = 0, how='any'), vmin = 0, vmax = 0.1, method='single')
    if ax == 'None': 
        matrix_map.get_figure().savefig(os.path.join(out_path, out_fname))
        plt.close()
    else:
        ax.set_title(out_fname.replace('.pdf',''), fontsize=24)
        return ax    
    

# create a function that compare a sample set with a given hit set
def sample_vs_hit(sample_notna, hit_notna, pop_notna, variable, test_type='ks', tail = 0, metric_cutoff = None, p_cutoff = 0.05):
    if test_type == 'ks':
        test_stat, test_pvalue = stats.ks_2samp(sample_notna, pop_notna)
    elif test_type == 'ttest':
        test_stat, test_pvalue = stats.ttest_1samp(sample_notna, pop_notna)
    elif test_type == 'binary':
        hit_stat = sum(hit_notna >= metric_cutoff)
        test_stat = sum(sample_notna >= metric_cutoff)
        oddsratio, test_pvalue = stats.fisher_exact([[test_stat, len(sample_notna)-test_stat], [hit_stat, len(hit_notna) - hit_stat]])            
    out = [test_stat, test_pvalue, np.mean(sample_notna)]
    return out


# create a function that generate n_sample random samples
def random_samples(df_pop, sample_size, id_coln = 'N. furzeri Protein Id', weights=None, \
                   n_sample = 10000, test_type='ks', tail = 0, metric_cutoff = None, p_cutoff = 0.05):
    out = dict()
    # create n_samples with or without weights
    for i in range(n_sample):
        # will sample from population without replacement
        if p_weights is None:
            sample_i = np.random.choice(df_pop[id_coln], sample_size, replace=False)
        else:
            sample_i = np.random.choice(df_pop[id_coln].values, sample_size, p=weights, replace=False)
        out[i] = sample_i
    return pd.DataFrame.from_dict(out, orient='columns')

# create a function that do permutation test and compare the # of times the sample_mean is different from all means
# http://t-redactyl.io/blog/2015/10/two-group-hypothesis-testing-permutation-tests.html     
def permutation_weights(df_hit, df_pop, variables, tissue, metric, fname_out=None, fpath_out=None, id_coln = 'N. furzeri Protein Id', df_samples = None,\
                        p_weights=None, n_sample = 10000, test_type='ks', tail = 0, metric_cutoff = None, p_cutoff = 0.05):
    # if you want weight population, please feed the pre-weight one
    # create n_samples with or without weights
    if df_samples == None:
        df_samples = random_samples(df_pop, df_hit.shape[0], weights=p_weights, n_sample=n_sample, \
                                test_type=test_type, tail=tail, metric_cutoff=metric_cutoff, p_cutoff=p_cutoff)
        df_samples.to_csv(os.path.join(fpath_out, '%s_%s.csv'%(fname_out, tissue)), index=False)
    # test enrichment for specific properties 
    out_obj = dict()
    out_result = np.zeros((len(variables), 8))
    for var_i, variable in enumerate(variables):
        test_matrix = np.zeros((n_sample,3))  # column: test-stat, p-val, mean_diff_to_pop 
        test_matrix[:] = np.nan
        hit_var = df_hit[variable].values
        pop_var = df_pop[variable].values
        hit_notna = hit_var[~np.isnan(hit_var)]
        pop_notna = pop_var[~np.isnan(pop_var)]
        if variable == 'NLLR': metric_cutoff = 0
        if variable == 'Disordered Frac': metric_cutoff = 0.3  
        if test_type == 'ks':
            hit_stat, hit_pvalue = stats.ks_2samp(hit_notna, pop_notna)
        elif test_type == 'ttest':
            hit_stat, hit_pvalue = stats.ttest_1samp(hit_notna, pop_notna)
        elif test_type == 'binary': 
            hit_stat = sum(hit_notna>= metric_cutoff)
            pop_pos = sum(pop_notna>= metric_cutoff)
            oddsratio, hit_pvalue = stats.fisher_exact([[pop_pos, len(pop_notna)- pop_pos], [hit_stat, df_hit.shape[0] - hit_stat]])   
        for i, s_i in enumerate(df_samples.columns):
            sample_i = df_samples[s_i]
            sample_var = df_pop[variable][df_pop[id_coln].isin(sample_i)].values
            sample_notna = sample_var[~np.isnan(sample_var)]         
            out = sample_vs_hit(sample_notna, hit_notna, pop_notna, variable, test_type=test_type, \
                                tail=tail, metric_cutoff=metric_cutoff, p_cutoff=p_cutoff)
            test_matrix[i,:] = out
            if out[0] == np.nan or out[1] == np.nan:
                print ('np.nan') 
        if hit_stat > np.mean(test_matrix[:,0]):
            p = sum(test_matrix[:,0] >= hit_stat)*1./n_sample
        else:
            p = sum(test_matrix[:,0] <= hit_stat)*1./n_sample
        if tail != 0:
            sig = True if p<= p_cutoff else False
        else:
            sig = True if p<= p_cutoff/2. else False
        z_score = (np.mean(hit_notna)- np.mean(test_matrix[:,2]))/np.std(test_matrix[:,2],ddof=0)
        out_obj['%s_%s_%s'%(tissue, metric, variable)] = {'hit.stat':hit_stat, 'hit.pval':hit_pvalue, 'hit.mean.zscore':z_score,\
               'hit.mean':np.mean(hit_notna), 'hit.size':len(hit_notna), 'pop.size':len(pop_notna),\
               'permutation.pval':p, 'test.type':test_type, 'permutation.tail':tail,\
               'test.matrix':test_matrix, 'n.sample':n_sample, 'p.cutoff':p_cutoff, 'sig':sig,\
                'Tissue':tissue, 'metric':metric, 'hit':hit_type, 'null':null_type}
        out_result[var_i,:] = [hit_stat, hit_pvalue, np.mean(hit_notna), z_score, p, len(hit_notna), len(pop_notna), n_sample]
    #    print 'stat = %.4f, p = %.4f, mean_diff=%.4f, sig = %s, out = %.4f'%(hit_stat, hit_pvalue, np.mean(hit)-np.mean(pop), sig, out)
    return out_obj, out_result    


# function to obtain all raw MS data and properties for a given df
def add_metric(df_sample, df_metrics, key_coln = None):
#    if df_metrics == None: df_raw = pd.read_csv(os.path.join('Correlomics', 'AllDectected_AllMetrics.csv'))
#    df_metrics = df_raw[df_raw['Sample.Type'] == 'AGG']
    if key_coln == None: key_coln = ['Tissue', 'N. furzeri Protein Id', 'N. furzeri Final Symbol', 'Human']
    key_coln = [i for i in key_coln if i in df_sample.columns]
    key_coln = [i for i in key_coln if i in df_metrics.columns]
    df_sample_metrics = df_sample[key_coln].merge(df_metrics, on = key_coln, how = 'left')
    return df_sample_metrics


fname_biophysics = 'Killifish_DISOPRED_PLAAC.csv'
fname_kaks = 'Ka_Ks_values_Nfur-to-otherfish_with-XM-XP-Ids.csv'
fname_cider = 'Killifish_Cider.csv'
f_path = 'Correlomics'
hit_path = 'Correlomics/Comp'
if not os.path.exists(hit_path): os.makedirs(hit_path)
tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis']

# # ######### Section 1: Combine all the different properties into one table for future import ########
# ## part 1. combine all the properties if it does not exist
# ## import the protein property files
# fname_biophysics = 'Killifish_DISOPRED_PLAAC.csv'
# # evolutionary results from param
# fname_kaks = 'Ka_Ks_values_Nfur-to-otherfish_with-XM-XP-Ids.csv'
# # localCIDER results
# fname_cider = 'Killifish_Cider.csv'
# fpath = 'Properties'
# if os.path.exists(os.path.join(fpath,'Killifish_DISOLLRCider.csv')) == False:
#     df_biophysics = pd.read_csv(os.path.join(fpath, fname_biophysics))
#     df_biophysics.rename(index=str, columns = {'Protein Id':'N. furzeri Protein Id'}, inplace = True)
#     df_kaks = pd.read_csv(os.path.join(fpath, fname_kaks))
#     df_cider = pd.read_csv(os.path.join(fpath, fname_cider))
#     df_cider.rename(index=str, columns = {'Protein Id':'N. furzeri Protein Id','Human':'Name'}, inplace = True)
#     df_kaks.rename(index=str, columns = {'XP id':'N. furzeri Protein Id'}, inplace = True)
#     df_metrics = pd.merge(df_biophysics, df_kaks.iloc[:,2:], how = 'outer', left_on = 'N. furzeri Protein Id', right_on = 'N. furzeri Protein Id')
#     df_metrics = pd.merge(df_metrics, df_cider.iloc[:,1:], how = 'outer', left_on = ['N. furzeri Protein Id','Seq'], right_on = ['N. furzeri Protein Id','Seq'])
#     df_metrics['Frac_QN'] = df_metrics['Frac_Q'] + df_metrics['Frac_N']
#     df_metrics['Frac_aliphatic'] = df_metrics['Frac_A'] + df_metrics['Frac_V']+ df_metrics['Frac_L'] + df_metrics['Frac_I']+ df_metrics['Frac_M']
#     df_metrics['Frac_aromatic'] = df_metrics['Frac_F'] + df_metrics['Frac_Y']+ df_metrics['Frac_W']
#     df_metrics['Frac_neg'] = df_metrics['Frac_D'] + df_metrics['Frac_E']
#     df_metrics['Frac_pos'] = df_metrics['Frac_H'] + df_metrics['Frac_K']+ df_metrics['Frac_R']
#     df_metrics['Frac_neu'] = df_metrics['Frac_S'] + df_metrics['Frac_T']+ df_metrics['Frac_N']+ df_metrics['Frac_Q']
#     df_metrics['Frac_FInumaa'] = df_metrics['FInumaa'] / df_metrics['AA Length']
#     df_metrics['Frac_FImaxrun'] = df_metrics['FImaxrun'] / df_metrics['AA Length']
#     df_metrics.to_csv(os.path.join(fpath, 'Killifish_DISOLLRCider.csv'), index = False)
# else:
#     df_metrics = pd.read_csv(os.path.join(fpath, 'Killifish_DISOLLRCider.csv'))
# ###  part 2. combine protein sequence properties with MS output
# fpath = 'Properties'
# if os.path.exists(os.path.join(fpath,'kf_combined_prop_metrics.csv')) == False:
#     df_metrics = pd.read_csv(os.path.join('Properties','Killifish_DISOLLRCider.csv'))
#     df_hit = pd.read_csv(os.path.join('', 'kf_combined_prop.csv'))
#     df_data = df_hit.merge(df_metrics, on = 'N. furzeri Protein Id', how='left')
#     df_data.to_csv(os.path.join('Properties','kf_combined_prop_metrics.csv'),index=False)
# else:
#     df_data = pd.read_csv(os.path.join('Properties','kf_combined_prop_metrics.csv'))
# ### part 2.1 compute a few metrics
# frac_NLLR = df_metrics['N. furzeri Protein Id'][df_metrics['NLLR']>=0].count()/df_metrics['N. furzeri Protein Id'].count()
# frac_diso = df_metrics['N. furzeri Protein Id'][df_metrics['Disordered Frac']>=0.3].count()/df_metrics['N. furzeri Protein Id'].count()
# df_nodup=df_data[['N. furzeri Protein Id','NLLR','Disordered Frac']].drop_duplicates()
# frac_MS_NLLR = df_nodup['N. furzeri Protein Id'][df_nodup['NLLR']>=0].count()/df_nodup['N. furzeri Protein Id'].count()
# frac_MS_diso = df_nodup['N. furzeri Protein Id'][df_nodup['Disordered Frac']>=0.3].count()/df_nodup['N. furzeri Protein Id'].count()
# ## part 2.2 plot delta vs. log2_Young_TL
# out_path = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join('Properties','AllDetected_metrics.csv'))
# tissues = ['Brain','Liver','Gut','Heart','Muscle','Skin','Testis']
# tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
# x="log2_Young_TL"; y="delta"
# fig,ax = plt.subplots(1,len(tissues),figsize=(len(tissues)*4,4), sharex=True, sharey=True)
# for i, tissue in enumerate(tissues):
#     df_data[(df_data['Sample.Type']=='AGG')&(df_data['Tissue']==tissue)].plot.scatter(x=x,y=y,\
#       color=tissues_colors[i],ax=ax[i], alpha=0.3)
# fig.savefig(os.path.join(out_path,'All_deltavlog2YTL.pdf'))
# hist_hue(df_data[df_data['Sample.Type']=='AGG'], ['delta'], 'AllDetected_deltahist', 'Tissue',\
#           hue_order = tissues,alpha = 1, bins = 50, normsum = True, directory = out_path)
# ######################################### End of Section 1 ########################################


# #####################################################################################################
# ##### Section 2:  Generate histogram of the all biophysical properties of a tissue
# ### update on 2020.04.23
# # initial input
# tissues_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
# tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
# panel_order = ['Frac_pos','FracNeg','FracCharged','FracExpanding','delta',\
#               'kappa','FracDisoPromoting', 'Disordered Frac', 'AA Length','FImaxrun',\
#               'PAPAprop','NLLR','Frac_neu','Frac_QN','MWScore',\
#               'MeanNetCharge','mean Ka/Ks','mean Ka','mean ks','pI',\
#               'Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic','log2_Young',\
#               'log2_Young_TL','OvY_logFC', 'OvY_prop_logFC']   # this is for sup figure
# hue_coln = 'Tissue'
# sample_type = 'AGG'#'TL'
# # read the killifish data
# fpath_kf = 'Correlomics'
# df_kf = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
# df_plot = df_kf[df_kf['Sample.Type']==sample_type]
# ### part 1: consider only the presence of a protein
# hist_coln_num = 5;
# hist_row_num = int(math.ceil(len(panel_order)/1./hist_coln_num))
# f, axes = plt.subplots(nrows = hist_row_num, ncols = hist_coln_num, figsize = (20, 20))
# for metric, ax in zip(panel_order, axes.flat):
#     ax = hist_samebin(metric, df_plot, ax, hue_coln, hue_order = tissues_order, colors_list=tissues_colors, bins= 50, normsum = True)
#     ax.set_title(metric)
# plt.tight_layout()
# f.savefig(os.path.join(fpath_kf,'%sProteomeHist.pdf'%sample_type))
# ### part 2: use the abundance of a protein as weight, then re-make the distribution, now do you see tissue-specific behavior?
# weight_coln = 'log2_Young_TL'#'Young_count_AGG'#'log2_Young'#'OvY_prop_logFC'#'log2_Young_TL'
# if weight_coln == 'Young_count_AGG': # create the prop weight_coln using young AGG data
#     df_plot.loc[:, weight_coln] = df_plot[['Young-1', 'Young-2', 'Young-3']].mean(axis=1)
# hist_coln_num = 5;
# hist_row_num = int(math.ceil(len(panel_order)/1./hist_coln_num))
# f, axes = plt.subplots(nrows = hist_row_num, ncols = hist_coln_num, figsize = (20, 20))
# for metric, ax in zip(panel_order, axes.flat):
#     ax = hist_samebin(metric, df_plot, ax, hue_coln, hue_order = tissues_order, colors_list=tissues_colors, bins= 50, normsum = True, weight_coln=weight_coln)
#     ax.set_title(metric)
# plt.tight_layout()
# f.savefig(os.path.join(fpath_kf,'%sProteomeHist_Weights_%s.pdf'%(sample_type, weight_coln)))
# #### part 3: plot cumulative plots
# weight_coln = 'Young_count_AGG'#'log2_Young_TL'#'OvY_prop_logFC'#'log2_Young_TL'
# if weight_coln == 'Young_count_AGG': # create the prop weight_coln using young AGG data
#     df_plot.loc[:, weight_coln] = df_plot[['Young-1', 'Young-2', 'Young-3']].mean(axis=1)
# hist_coln_num = 5;
# hist_row_num = int(math.ceil(len(panel_order)/1./hist_coln_num))
# f, axes = plt.subplots(nrows = hist_row_num, ncols = hist_coln_num, figsize = (20, 20))
# for metric, ax in zip(panel_order, axes.flat):
#     ax = hist_cum(metric, df_plot, ax, hue_coln, hue_order = tissues_order, colors_list=tissues_colors, bins= 50, normsum = True, weight_coln=weight_coln)
#     ax.set_title(metric)
# plt.tight_layout()
# f.savefig(os.path.join(fpath_kf,'%sProteomeCumHist_Weights_%s.pdf'%(sample_type, weight_coln)))
# ############################### End of Section 10 ##############################

    
# ###################################################################################################
# ###### Section 3: Get the correlation between evolutionary rates, expression data and disorder metrics etc. #####
######### 2020.04.20 plot ########
# this aim to compute the spearman r of hit (z-score filtered or not) vs. null background and show the results in heatmap format
# tissues_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
# fpath_all = 'YCMSResults'
# fpath_out = 'Correlomics/Evolution'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# ## read the input file
# df_all = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
# # input the age comparison of interest
# x_list = ['log2_Young_TL', 'log2_Young', 'log2_Young_prop','OvY_logFC', 'OvY_prop_logFC']
# y_list = ['mean ks', 'mean Ka', 'mean Ka/Ks']
# metric_list = ['All','AGG','Prop']
# p_cutoff=0.05
# z_cutoff=None
# age1 = 'Young'; age2 = 'Old'
# comp = '%sv%s'%(age2[0], age1[0])
# tissues_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin']
# if age2 != 'Tert': tissues_order.append('Testis')
# ###### part 1: create the matrix to store the spearman correlation rho
# row_pairs = [(x_i, y_i) for y_i in y_list for x_i in x_list]
# row_list = ['%s.VS.%s'%(y_i, x_i) for y_i in y_list for x_i in x_list]
# coln_list = []
# for metric in metric_list:
#     for i in (tissues_order+['all']):
#         if metric == 'All':
#             coln_list.append('DetectedAGG_%s'%i)
#         else:
#             if z_cutoff == None:
#                 coln_list.append('%sSigPos_%s_%s'%(metric,comp,i))
#             else:
#                 coln_list.append('%sSigPosZ%s_%s_%s'%(metric,z_cutoff,comp,i))
# rho_matrix = np.empty((len(row_list),len(coln_list)))
# rho_matrix[:] = np.nan
# pval_matrix = np.empty((len(row_list),len(coln_list)) )
# pval_matrix[:] = np.nan
# for m_c, metric in enumerate(metric_list):
#     if metric  == 'Prop':
#         fc_coln = '%sv%s_prop_logFC'%(age2[0], age1[0])
#         pval_coln = '%sv%s_prop_pval'%(age2[0], age1[0])
#     else:
#         fc_coln = '%sv%s_logFC'%(age2[0], age1[0])
#         pval_coln = '%sv%s_pval'%(age2[0], age1[0])
#     z_coln = '%s_zscore'%fc_coln
#     df_test = df_all[df_all['Sample.Type']=='AGG']
#     if metric !='All':
#         df_test = df_test[(df_test[pval_coln]<=p_cutoff) & (df_test[fc_coln]>=0)] if z_cutoff == None else df_test[(df_test[pval_coln]<=p_cutoff) & (df_test[z_coln]>=z_cutoff)]
#     for i, pair in enumerate(row_pairs):
#         x_i, y_i = pair
#         for t_c, tissue in enumerate(tissues_order):
#             coln_index = m_c*(len(tissues_order)+1) + t_c
#             rho, pval = scipy.stats.spearmanr(df_test[[x_i,y_i]][df_test['Tissue']==tissue].dropna())
#             rho_matrix[i, coln_index] = rho
#             pval_matrix[i, coln_index] = pval
#         rho_all, pval_all = scipy.stats.spearmanr(df_test[[x_i,y_i]].dropna())
#         rho_matrix[i, coln_index+1] = rho_all
#         pval_matrix[i, coln_index+1] = pval_all
# df_rho = pd.DataFrame(data = rho_matrix, index = row_list, columns = coln_list)
# df_pval = pd.DataFrame(data = pval_matrix, index = row_list, columns = coln_list)
# df_rho.to_csv(os.path.join(fpath_out,'SigPos_%s_spearmanrho.csv'%comp)) if z_cutoff == None else df_rho.to_csv(os.path.join(fpath_out,'SigPosZ%s_%s_spearmanrho.csv'%(z_cutoff,comp)))
# df_pval.to_csv(os.path.join(fpath_out,'SigPos_%s_spearmanpval.csv'%comp)) if z_cutoff == None else df_pval.to_csv(os.path.join(fpath_out,'SigPosZ%s_%s_spearmanpval.csv'%(z_cutoff,comp)))
# ###### part 2: visualize some scatter plots
# for metric in metric_list:
#     df_test = df_all[df_all['Sample.Type']=='AGG']
#     if metric  == 'Prop':
#         fc_coln = '%sv%s_prop_logFC'%(age2[0], age1[0])
#         pval_coln = '%sv%s_prop_pval'%(age2[0], age1[0])
#     else:
#         fc_coln = '%sv%s_logFC'%(age2[0], age1[0])
#         pval_coln = '%sv%s_pval'%(age2[0], age1[0])
#     z_coln = '%s_zscore'%fc_coln
#     if metric !='All':
#         df_test = df_test[(df_test[pval_coln]<=p_cutoff) & (df_test[fc_coln]>=0)] if z_cutoff == None else df_test[(df_test[pval_coln]<=p_cutoff) & (df_test[z_coln]>=z_cutoff)]
#     fig_row = len(x_list)
#     fig_coln = len(tissues_order)+1
#     f, axes = plt.subplots(nrows = fig_row, ncols = fig_coln, figsize = (fig_coln*4, fig_row*4), sharey=True)
#     for i, x in enumerate(x_list):
#         for j, tissue in enumerate(tissues_order):
#             df_plot = df_test[df_test['Tissue']==tissue]
#             sns.regplot(x = x, y ='mean Ka', data = df_plot, ax = axes[i,j], color='k',scatter_kws={'alpha':0.2},line_kws={"color": "red"})
#             rho, pval = scipy.stats.spearmanr(df_plot[[x,'mean Ka']].dropna())
#             axes[i,j].xaxis.label.set_size(16)
#             axes[i,j].yaxis.label.set_size(16)
#             axes[i,j].tick_params(labelsize=12)
#             axes[i,j].annotate('%s r = %.2f, p = %.2f'%(tissue, rho, pval), xy = (0.05, 0.9), xycoords='axes fraction',fontsize=14)
#         df_plot = df_all[df_all['Sample.Type']=='AGG']
#         sns.regplot(x = x, y ='mean Ka', data = df_plot, ax = axes[i,j+1], color='k',scatter_kws={'alpha':0.2}, line_kws={"color": "red"})
#         rho, pval = scipy.stats.spearmanr(df_plot[[x,'mean Ka']].dropna())
#         axes[i,j+1].xaxis.label.set_size(16)
#         axes[i,j+1].yaxis.label.set_size(16)
#         axes[i,j+1].tick_params(labelsize=12)
#         axes[i,j+1].annotate('All Tissues r = %.2f, p = %.2f'%(rho, pval), xy = (0.05, 0.9), xycoords='axes fraction',fontsize=14)
#     if metric == 'All':
#         f.savefig(os.path.join(fpath_out,'%s_meanKa_scatter.png'%metric))
#     else:
#         f.savefig(os.path.join(fpath_out,'%sSigPos_%s_meanKa_scatter.png'%(metric,comp))) if z_cutoff == None else f.savefig(os.path.join(fpath_out,'%sSigPosZ%s_%s_meanKa_scatter.png'%(metric,z_cutoff,comp)))
#     plt.close(f)
# ###### part 3: is brain aggregate unique and its diso related to mean Ka and age-associated aggregation?
# fpath_out = 'Correlomics/Evolution/Pairwise'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# metric = 'All' # choose between 'All' and 'Prop'
# constant_y = 'NLLR'
# matrix_coln = ['Disordered Frac', 'FracDisoPromoting', 'NLLR', 'AA Length', 'mean Ka','delta','kappa','uversky_hydropathy','log2_Young_TL','OvY_prop_logFC', 'OvY_logFC']
# if constant_y  in matrix_coln: matrix_coln.remove(constant_y)
# fig_row = len(tissues_order)+1
# fig_coln = len(matrix_coln)
# f, axes = plt.subplots(nrows = fig_row, ncols = fig_coln, figsize = (fig_coln*4, fig_row*4), sharey=True) #do not sharex because the plot will shrink into dots
# df_plot = df_all[df_all['Sample.Type']=='AGG']
# for i, tissue in enumerate(tissues_order+['all']):
#     df_ij = df_plot[df_plot['Tissue']==tissue] if tissue !='all' else df_plot
#     if metric == 'Prop':
#         df_ij = df_ij[(df_ij['OvY_prop_logFC']>=0)&(df_ij['OvY_prop_pval']<=0.05)]
#     elif metric == 'AGG':
#         df_ij = df_ij[(df_ij['OvY_logFC']>=0)&(df_ij['OvY_pval']<=0.05)]
#     for j, variable in enumerate(matrix_coln):
#         sns.regplot(x = variable, y = constant_y, data = df_ij, ax = axes[i,j], color = 'k', scatter_kws={'alpha':0.2}, line_kws={"color": "red"})
#         rho, pval = scipy.stats.spearmanr(df_ij[[variable, constant_y]].dropna())
#         axes[i,j].xaxis.label.set_size(16)
#         axes[i,j].yaxis.label.set_size(16)
#         axes[i,j].tick_params(labelsize=12)
#         axes[i,j].annotate('%s r = %.2f, p = %.2f'%(tissue, rho, pval), xy = (0.05, 0.9), xycoords='axes fraction',fontsize=14)
# f.savefig(os.path.join(fpath_out,'%sSigPos_OvY_%sVSmetrics.pdf'%(metric, constant_y.replace(' ','')))) if metric !='All' else f.savefig(os.path.join(fpath_out,'%s_%sVSmetrics.pdf'%(metric, constant_y.replace(' ',''))))
# plt.close(f)
## comment: do not use grid_correlomics for this because Ka has lots of missing values so the other pairwise comparison will be very off!!!
# ##### End of section 3 #####
# ####################################################################################################
        

# ############ Section 4: Implement permutation test to compare hits (increase with age) vs. null background or other tissue hits
# input the age comparison of interest
# updated on 2020.06.30
# either use weight_coln = None (no weights) or weight_coln = 'Young_count' (with weights for sampling)
weight_coln = None ##None  #'Young_count'
query_metrics = ['Prop', 'TL'] if weight_coln is None else ['AGG'] #'AGG', 
age1 = 'Young'; age2 = 'Old'; comp = '%sv%s'%(age2[0], age1[0])
test_type = 'ks'#'binary'#'ks'
hit_sign = 'Pos' #'Neg'
tissues_order = ['Brain','Gut','Heart','Liver','Muscle','Skin'] if 'T' in comp else ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']
# specify the parameter you want to test
panel_order = ['FracPos', 'Frac_pos', 'FracNeg', 'FracCharged', 'FracExpanding', 'delta', 'deltaMax', 'kappa',\
                      'FracDisoPromoting', 'Disordered Frac', 'AA Length', 'Longest Disorder Stretch','Frac_FInumaa','Frac_FImaxrun',\
                      'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu', 'Frac_QN', 'MWScore', 'VITmaxrun', 'PRDscore', 'MeanNetCharge',\
                      'mean Ka/Ks','mean Ka','mean ks',\
                      'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
                      'log2_Young_TL','log2_Old_TL','log2_Young_prop','OvY_logFC', 'OvY_prop_logFC'] ### do not change the order here, i grouped it based on the pairwise-correlation and clustering
# panel_order = ['NLLR', 'Disordered Frac']
csv_colns = ['hit.stat','hit.pval','hit.mean','hit.mean.zscore','permutation.pval','hit.size','pop.size','n.sample']
fpath_in = 'Correlomics'
df_all = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
if weight_coln == 'Young_count':
    # create the prop weight_coln
    df_all.loc[:,weight_coln] = df_all[['Young-1', 'Young-2', 'Young-3']].mean(axis=1)
####### part 1.1 conduct permuation test to estimate z-score and p-value (with or without weights)
for metric in query_metrics:
    null_type = '%sAll_%s'%(metric,comp)
    df_null = pd.read_csv(os.path.join('%s/Comp'%fpath_in, '%s.csv'%null_type))
    if metric !='TL':
        df_null_metrics = df_all[df_all['Sample.Type']=='AGG']
    else:
        df_null_metrics = df_all[df_all['Sample.Type']==metric]
    p_cutoff=0.05
    for z_cutoff in [0.75, None]:
        fpath_out = '%s/%sHits_%s/PermutationWeights'%(fpath_in,metric,comp) if weight_coln is not None else '%s/%sHits_%s/Permutation'%(fpath_in,metric,comp)
        if not os.path.exists(fpath_out): os.makedirs(fpath_out)
        hit_type = '%sSig%s_%s'%(metric, hit_sign, comp) if z_cutoff is None else '%sSig%sZ%s_%s'%(metric, hit_sign, z_cutoff, comp)
        print (hit_type)
        fname_hit = '%sSig_%s.csv'%(metric,comp) if z_cutoff is None else '%sSigZ%s_%s.csv'%(metric,z_cutoff,comp)
        df_hit = pd.read_csv(os.path.join('%s/Comp'%fpath_in, fname_hit))
        # specify the parameter you want to test
        df_hit_metrics = add_metric(df_hit[df_hit['Hit Type']==hit_type], df_metrics=df_null_metrics)
        # permutation_obj = dict()
        permutation_result = np.zeros((len(panel_order)*len(tissues_order), len(csv_colns))) # first coln to test-stats, second to store test-pval, third to store hit_mean, fourth to store store z-score, fifth to store permutation-pval
        permutation_result[:,:] = np.nan
        for j, tissue in enumerate(tissues_order):
            df_pop = df_null_metrics[df_null_metrics['Tissue']==tissue]
            if weight_coln != None:
                df_weights = df_pop[weight_coln]
                df_weights = df_weights.fillna(df_weights.min())
                p_weights = df_weights.div(df_weights.sum())
            else:
                p_weights = None
            permutation_obj, permutation_result_j = permutation_weights(df_hit_metrics[df_hit_metrics['Tissue']==tissue], df_pop, panel_order,\
                                                    tissue, metric, fname_out=hit_type, fpath_out=fpath_out, id_coln = 'N. furzeri Protein Id',\
                                                    df_samples = None, p_weights = p_weights, n_sample = 10000,\
                                                    test_type = test_type, tail = 0, metric_cutoff = None, p_cutoff = 0.05)
            permutation_result[j*len(panel_order):(j+1)*len(panel_order),:] = permutation_result_j
        df_permutation = pd.DataFrame(data = permutation_result, columns = csv_colns)
        df_permutation['Tissue'] = [tissue for tissue in tissues_order for i in range(len(panel_order))]
        df_permutation['Metric'] = panel_order * len(tissues_order)
        df_permutation['test.type'] = '%s_to_pop'%test_type
        fout_prefix = '%s_%s_%s_weights'%(hit_type, test_type, weight_coln) if weight_coln is not None else '%s_%s'%(hit_type, test_type)
        df_permutation.to_csv(os.path.join(fpath_out, '%s.csv'%fout_prefix),index=False)
        f_obj = open(os.path.join(fpath_out, '%s.obj'%fout_prefix), 'wb')
        pickle.dump(permutation_obj, f_obj)
# #### part 2.1 load specific f_obj of interest from the permuation run to look at distribution
# metric = 'Prop'
# comp = 'OvY'
# z_cutoff = None # 0.75
# hit_sign = 'Pos'
# hit_type = '%sSig%s_%s'%(metric, hit_sign, comp)  if z_cutoff == None else '%sSig%sZ%s_%s'%(metric, hit_sign, z_cutoff, comp)
# f_prefix = '%s_ks_YAGGweights'%hit_type #'%s_%s'%(hit_type, test_type) for the ones without weights
# fpath_out = 'Correlomics/%sHits_%s/Permutation'%(metric, comp)
# f_obj = open(os.path.join(fpath_out, '%s_ks.obj'%hit_type),'rb')
# permutation_obj = pickle.load(f_obj)
############################### End of Section 4 ##############################
    

###### Section 5: Test normality of metric distribution, then that with added abundance weights 
##input the age comparison of interest
# test_type = 'ks'
# tissues_order = ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']
# panel_order = ['FracPos','Frac_pos','FracNeg','FracCharged','FracExpanding','delta','deltaMax','kappa',\
#                       'FracDisoPromoting', 'Disordered Frac', 'AA Length',\
#                       'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu','Frac_QN','MWScore','VITmaxrun', 'PRDscore','MeanNetCharge',\
#                       'mean Ka/Ks','mean Ka','mean ks',\
#                       'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
#                       'log2_Young_TL','log2_Old_TL','log2_Young_prop','OvY_logFC', 'OvY_prop_FC'] ### do not change the order here, i grouped it based on the pairwise-correlation and clustering
# fpath_in = 'Correlomics'
# fpath_out = 'Correlomics/Model/Normality'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# test_type = 'ks'#'anderson'#'shapiro' #'ks'
# p_cutoff = 0.05
# df_all = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
# #### part 1: do this with no weights assigned to the values 
# for sample_type in ['AGG','TL']:
#     norm_tests = np.zeros((len(panel_order)*len(tissues_order),3))
#     norm_tests[:,:] = np.nan 
#     for j, tissue in enumerate(tissues_order):
#         for i, variable in enumerate(panel_order):
#             df_test = df_all[variable][(df_all['Tissue']==tissue) & (df_all['Sample.Type']==sample_type)].dropna()
#             if test_type == 'ks':
#                 norm_tests[i+j*len(panel_order),0:2] = stats.kstest(df_test.values, 'norm')
#             elif test_type == 'shapiro':
#                 norm_tests[i+j*len(panel_order),0:2] = stats.shapiro(df_test.values)
#             elif test_type == 'anderson':
#                 anderson_out = stats.anderson(df_test.values) #normal/exponenential, 15%, 10%, 5%, 2.5%, 1%
#                 norm_tests[i+j*len(panel_order),0:2] = anderson_out[0], anderson_out[1][2]
#             norm_tests[i+j*len(panel_order),2] = df_test.shape[0]
#     if test_type !='anderson':
#         df_norm = pd.DataFrame(data = norm_tests, columns = ['%s.stat'%test_type, '%s.pval'%test_type, 'n']) 
#     else:
#         df_norm = pd.DataFrame(data = norm_tests, columns = ['%s.stat'%test_type, '%s_p%s.stat'%(test_type,p_cutoff), 'n']) 
#     df_norm['Tissue'] = [tissue for tissue in tissues_order for i in range(len(panel_order))]
#     df_norm['Metric'] = panel_order * len(tissues_order)
#     df_norm.to_csv(os.path.join(fpath_out, '%s_norm_%s.csv'%(sample_type,test_type)),index=False)
# ### part 2: do this with weights from log2_Young_TL or log2_Young_AGG assigned to the values 
# weight_coln = 'log2_Young'#'OvY_prop_logFC'#'log2_Young_TL'
# for sample_type in ['AGG','TL']:
#     norm_tests = np.zeros((len(panel_order)*len(tissues_order),3))
#     norm_tests[:,:] = np.nan 
#     for j, tissue in enumerate(tissues_order):
#         for i, variable in enumerate(panel_order):
#             if variable == weight_coln:
#                 continue
#             df_ij = df_all[[variable,weight_coln]][(df_all['Tissue']==tissue) & (df_all['Sample.Type']==sample_type)]
#             df_test = df_ij[variable]*df_ij[weight_coln]
#             df_test = df_test.dropna()
#             if test_type == 'ks':
#                 norm_tests[i+j*len(panel_order),0:2] = stats.kstest(df_test.values, 'norm')
#             elif test_type == 'shapiro':
#                 norm_tests[i+j*len(panel_order),0:2] = stats.shapiro(df_test.values)
#             elif test_type == 'anderson':
#                 anderson_out = stats.anderson(df_test.values) #normal/exponenential, 15%, 10%, 5%, 2.5%, 1%
#                 norm_tests[i+j*len(panel_order),0:2] = anderson_out[0], anderson_out[1][2]
#             norm_tests[i+j*len(panel_order),2] = df_test.shape[0]
#     if test_type !='anderson':
#         df_norm = pd.DataFrame(data = norm_tests, columns = ['%s.stat'%test_type, '%s.pval'%test_type, 'n']) 
#     else:
#         df_norm = pd.DataFrame(data = norm_tests, columns = ['%s.stat'%test_type, '%s_p%s.stat'%(test_type,p_cutoff), 'n']) 
#     df_norm['Tissue'] = [tissue for tissue in tissues_order for i in range(len(panel_order))]
#     df_norm['Metric'] = panel_order * len(tissues_order)
#     df_norm.to_csv(os.path.join(fpath_out, '%s_norm_%s_YAGGweights.csv'%(sample_type,test_type)),index=False)  
################################ End of Section 5 ##############################
    

# # ######### Section 6: feature of proteins that are aggregation-prone in young ######
# weight_coln = 'count_Young_TL' #'count_Old_TL'#None #'count_Young_TL'
# query_metrics = ['AGG','Prop'] if weight_coln is None else ['AGG']
# age = 'Young'
# test_type = 'ks'#'binary'#'ks' # use binary for panel_order = ['NLLR', 'Disordered Frac'], use ks for the longer list
# hit_sign = 'Pos' #'Neg'
# z_cutoff = 2
# direction = 'greater'
# tissues_order = ['Brain','Gut','Heart','Liver','Muscle','Skin'] if 'T' in age else ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']
# # specify the parameter you want to test
# panel_order = ['FracPos','Frac_pos','FracNeg','FracCharged','FracExpanding','delta','deltaMax','kappa',\
#                       'FracDisoPromoting', 'Disordered Frac','Longest Disorder Stretch', 'AA Length', 'Frac_FInumaa','Frac_FImaxrun',\
#                       'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu','Frac_QN','MWScore','VITmaxrun', 'PRDscore','MeanNetCharge',\
#                       'mean Ka/Ks','mean Ka','mean ks',\
#                       'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
#                       'log2_Young_TL','log2_Old_TL','log2_Young_prop','OvY_logFC', 'OvY_prop_logFC'] ## use this if metric is TL####
# # panel_order = ['NLLR', 'Disordered Frac']
# csv_colns = ['hit.stat','hit.pval','hit.mean','hit.mean.zscore','permutation.pval','hit.size','pop.size','n.sample']
# fpath_in = 'Correlomics'
# df_all = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
# ######## part 1.1 conduct permuation test to estimate z-score and p-value (with or without weights)
# for metric in query_metrics:#['TL', 'AGG', 'Prop']:
#     null_type = '%sAll_OvY'%metric
#     df_null = pd.read_csv(os.path.join('%s/Comp'%fpath_in, '%s.csv'%null_type))
#     if metric !='TL':
#         df_null_metrics = df_all[df_all['Sample.Type']=='AGG']
#     else:
#         df_null_metrics = df_all[df_all['Sample.Type']==metric]
#     hit_path = 'Correlomics/AgeSpecific'
#     hit_fname = '%s_%sZ%s'%(age, metric, z_cutoff)
#     hit_type = '%sPosZ%s'%(metric, z_cutoff) if z_cutoff >0 else '%sNegZ%s'%(metric, z_cutoff)
#     df_hit = pd.read_csv(os.path.join(hit_path, '%s.csv'%hit_fname))
#     out_path = '%s/Permutation'%hit_path if weight_coln == None else '%s/PermutationWeights'%hit_path
#     if not os.path.exists(out_path): os.makedirs(out_path)
#     # specify the parameter you want to test
#     df_hit_metrics = add_metric(df_hit[df_hit['Hit Type']==hit_type], df_metrics=df_null_metrics)
#     # permutation_obj = dict()
#     permutation_result = np.zeros((len(panel_order)*len(tissues_order),len(csv_colns))) # first coln to test-stats, second to store test-pval, third to store hit_mean, fourth to store store z-score, fifth to store permutation-pval
#     permutation_result[:,:] = np.nan
#     for j, tissue in enumerate(tissues_order):
#         df_pop = df_null_metrics[df_null_metrics['Tissue']==tissue]
#         if weight_coln != None:
#             df_weights = df_pop[weight_coln]
#             df_weights = df_weights.fillna(df_weights.min())
#             p_weights = df_weights.div(df_weights.sum())
#         else:
#             p_weights = None
#         permutation_obj, permutation_result_j = permutation_weights(df_hit_metrics[df_hit_metrics['Tissue']==tissue], df_pop, panel_order,\
#                                                 tissue, metric, fname_out=hit_type, fpath_out=out_path, id_coln = 'N. furzeri Protein Id',\
#                                                 df_samples = None, p_weights = p_weights, n_sample = 10000, test_type = test_type)
#         permutation_result[j*len(panel_order):(j+1)*len(panel_order),:] = permutation_result_j
#     df_permutation = pd.DataFrame(data = permutation_result, columns = csv_colns)
#     df_permutation['Tissue'] = [tissue for tissue in tissues_order for i in range(len(panel_order))]
#     df_permutation['Metric'] = panel_order * len(tissues_order)
#     df_permutation['test.type'] = '%s_to_pop'%test_type
#     out_fname = '%s_%s'%(hit_fname, test_type) if weight_coln is None else '%s_%s_%s_weights'%(hit_fname, test_type, weight_coln)
#     df_permutation.to_csv(os.path.join(out_path, '%s.csv'%out_fname),index=False)
#     f_obj = open(os.path.join(out_path, '%s.obj'%out_fname), 'wb')
#     pickle.dump(permutation_obj, f_obj)
#### part 2.1 load specific f_obj of interest from the permuation run to look at distribution
# p_cutoff = .05
# z_cutoff = 2
# age = 'Young'
# weights = 'count_Young_TL' #'count_Young_TL'#None #'count_Old_TL'
# metric ='AGG'
# # load the object
# fpath_hit = 'Correlomics/AgeSpecific/Permutation' if weights == None else 'Correlomics/AgeSpecific/PermutationWeights' 
# fname_hit = '%s_%sZ%s_ks'%(age, metric, z_cutoff) if weights == None else '%s_%sZ%s_ks_%sweights'%(age, metric, z_cutoff, weights) # '%sSig%s_%s_ks'%(metric, hit_dir, age_comp)
# f_obj = open(os.path.join(fpath_hit, '%s.obj'%fname_hit),'rb')
# permutation_obj = pickle.load(f_obj)
# ##use the following lines to look at individual parameters
# permutation_obj['Brain_AGG_delta']
################################# End of Section 6 ##############################   
    

#######################################################################################################
########## Section 7: for supp figure S3 (tweak) #####
############ 2020.05.03 plot ########
## this aim to compute the spearman r of AGG and Prop (z-score filtered or not) from a particular age vs. different features
# # input the age comparison of interest
# age = 'TERT'
# ## read the input file
# df_data = pd.read_csv(os.path.join('Properties', 'kf_combined_prop_metrics.csv'))
# fpath_out = 'Correlomics/AgeSpecific/OneVariableRegression'
# if not os.path.exists(fpath_out): os.makedirs(fpath_out)
# # specify the parameter you want to test
# panel_order = ['FracPos','Frac_pos','FracNeg','FracCharged','FracExpanding','delta','deltaMax','kappa',\
#               'FracDisoPromoting', 'Disordered Frac', 'AA Length',\
#               'FInumaa','FImaxrun', 'PAPAprop','NLLR','Frac_neu','Frac_QN','MWScore','VITmaxrun', 'PRDscore','MeanNetCharge',\
#               'mean Ka/Ks','mean Ka','mean ks',\
#               'pI','Frac_aromatic','Omega','uversky_hydropathy','Frac_aliphatic',\
#               'log2_Young_TL','log2_Old_TL'] ### do not change the order here, i grouped it based on the pairwise-correlation and clustering
# tissues_order = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin'] if age == 'TERT' else ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
# for metric in ['AGG', 'Prop']:
#     p_cutoff=0.05
#     for z_cutoff in [None, 2]:
#         if metric == 'Prop':
#             value = 'log2_%s_prop'%age
#             df_all = df_data[df_data['Sample.Type']=='AGG']
#         else:
#             value = 'log2_%s'%age
#             df_all = df_data[df_data['Sample.Type']==metric]
#         z_coln = '%s_zscore'%value
#         # create the matrix to store the spearman correlation rho
#         rho_all_matrix = np.empty((len(panel_order),len(tissues_order)+1))
#         rho_pos_matrix = np.empty((len(panel_order),len(tissues_order)+1))
#         rho_all_matrix[:] = np.nan;
#         rho_pos_matrix[:] = np.nan;
#         pval_all_matrix = np.empty((len(panel_order),len(tissues_order)+1))
#         pval_pos_matrix = np.empty((len(panel_order),len(tissues_order)+1))
#         pval_all_matrix[:] = np.nan;
#         pval_pos_matrix[:] = np.nan;
#         for x_i, variable in enumerate(panel_order):
#             for y_i, tissue in enumerate(tissues_order):
#                 df_test = df_all[df_all['Tissue']==tissue]
#                 rho_all, pval_all = stats.spearmanr(df_test[[variable,value]].dropna())
#                 if z_cutoff == None:
#                     rho_pos, pval_pos = stats.spearmanr(df_test[[variable,value]][df_test[value]>=0].dropna())
#                 else:
#                     rho_pos, pval_pos = stats.spearmanr(df_test[[variable,value]][df_test[z_coln]>=z_cutoff].dropna())
#                 rho_all_matrix[x_i,y_i] = rho_all
#                 rho_pos_matrix[x_i,y_i] = rho_pos
#                 pval_all_matrix[x_i,y_i] = pval_all
#                 pval_pos_matrix[x_i,y_i] = pval_pos
#             df_test = df_all[df_all['Sample.Type']=='AGG']
#             rho_all, pval_all = stats.spearmanr(df_test[[variable,value]].dropna())
#             if z_cutoff == None:
#                 rho_pos, pval_pos = stats.spearmanr(df_test[[variable,value]][df_test[value]>=0].dropna())
#             else:
#                 rho_pos, pval_pos = stats.spearmanr(df_test[[variable,value]][df_test[value]>=z_cutoff].dropna())
#             rho_all_matrix[x_i,y_i+1] = rho_all
#             rho_pos_matrix[x_i,y_i+1] = rho_pos
#             pval_all_matrix[x_i,y_i+1] = pval_all
#             pval_pos_matrix[x_i,y_i+1] = pval_pos
#         df_rho_all = pd.DataFrame(data = rho_all_matrix, index = panel_order, columns = tissues_order+['all'])
#         df_rho_pos = pd.DataFrame(data = rho_pos_matrix, index = panel_order, columns = tissues_order+['all'])
#         df_pval_all = pd.DataFrame(data = pval_all_matrix, index = panel_order, columns = tissues_order+['all'])
#         df_pval_pos = pd.DataFrame(data = pval_pos_matrix, index = panel_order, columns = tissues_order+['all'])
#         df_rho_all.to_csv(os.path.join(fpath_out,'%sAll_%s_spearmanrho.csv'%(metric,age)))
#         df_rho_pos.to_csv(os.path.join(fpath_out,'%sPos_%s_spearmanrho.csv'%(metric,age))) if z_cutoff == None else df_rho_pos.to_csv(os.path.join(fpath_out,'%sPosZ%s_%s_spearmanrho.csv'%(metric,z_cutoff,age)))
#         df_pval_all.to_csv(os.path.join(fpath_out,'%sAll_%s_spearmanpval.csv'%(metric,age)))
#         df_pval_pos.to_csv(os.path.join(fpath_out,'%sPos_%s_spearmanpval.csv'%(metric,age))) if z_cutoff == None else df_pval_pos.to_csv(os.path.join(fpath_out,'%sPosZ%s_%s_spearmanpval.csv'%(metric,z_cutoff,age)))
# ################################## End of Section 7 ##############################
#####################################################################################################
