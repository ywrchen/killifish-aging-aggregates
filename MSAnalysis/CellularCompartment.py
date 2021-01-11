#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:00:55 2019

this script is used to analyze the quantify the cellular localization changes in aggregate
@author: yiwenchen
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math


mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams.update({'font.size': 40})
mpl.rcParams.update({'font.size': 8})

# plot pie chart for each compartment for a given tissue, each compartment is plot on a separate subplot
def pie_CM(df_CM, st, tissue, cutoff, variable, compartments, cutoff_type ='per', p_coln = None, p_cutoff = 0.05, rel = False):
    # here cutoff_type can either be based on percentile as before "per", or based on "value" of the variable of interest (i.e. "fc" or "z-score")
    # cutoff = [0,100] when cutoff_type == 'per', or else cutoff = [upper_cutoff, lower_cutoff]
    fig, axes = plt.subplots(nrows = 5, ncols = 4, figsize=(100, 120), sharex = True, sharey = True)
    ax_list = axes.flat
    df_tissue = df_CM[(df_CM['Tissue'] == tissue) & (df_CM['Sample.Type'] == st)]
    if cutoff_type == 'per':
        tissue_top = df_tissue[df_tissue[variable] >= df_tissue[variable].quantile(cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
        tissue_bottom = df_tissue[df_tissue[variable] <= df_tissue[variable].quantile(1-cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
    elif cutoff_type == 'value':
        if p_coln == None:
            tissue_top = df_tissue[df_tissue[variable] >= cutoff[0]].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            tissue_bottom = df_tissue[df_tissue[variable] <= cutoff[1]].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
        else:
            tissue_top = df_tissue[(df_tissue[variable]>=cutoff[0]) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            tissue_bottom = df_tissue[(df_tissue[variable]<= cutoff[1]) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
    tissue_all = df_tissue.groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
    colors = ['gold', 'lightskyblue', 'gray']
    for i, cm_i in enumerate(compartments):
        if (cm_i not in tissue_top.index):
            up = 0
        else:
            up = tissue_top.loc[cm_i].values[0]
        if (cm_i not in tissue_bottom.index):
            down = 0
        else:
            down = tissue_bottom.loc[cm_i].values[0]
        if up != 0 or down != 0:
            tot = tissue_all.loc[cm_i].values[0]
            data = [up, down, tot-up-down]
            if rel == True:
                wedges, texts, autotexts = ax_list[i].pie(data, colors = colors, radius = math.sqrt(tot)*0.07, autopct = lambda pct: '%.2f%%'%pct)
            else:
                wedges, texts, autotexts = ax_list[i].pie(data, colors = colors, autopct = lambda pct: '%.2f%%'%pct)
            if cutoff_type == 'per':
                ax_list[i].legend(wedges, ['>=%s percentile'%int(cutoff*100), '<= %s percentile'%(100-100*cutoff), 'other'], title = 'fold change', loc='center left')
            elif cutoff == 'value':
                ax_list[i].legend(wedges, ['>=%s'%cutoff[0], '<= %s'%cutoff[1], 'other'], title = '%s value'%variable, loc='center left')
            ax_list[i].set_title('%s %s (%.d)'%(tissue, cm_i, tot))
    return fig, axes


# do accounting on number of proteins that are in top or bottom 10 percentile in each tissue, then figure out where they fall
# in terms of cellular comparment
def tissue_CM(df_CM, st, variable, cutoff, cutoff_type='per', p_coln = None, p_cutoff =0.05):
    # here cutoff_type can either be based on percentile as before "per", or based on "value" of the variable of interest (i.e. "fc" or "z-score")
    # cutoff = [0,100] when cutoff_type == 'per', or else cutoff = [upper_cutoff, lower_cutoff]
    frames = []
    tissues = ['Brain', 'Gut', 'Heart', 'Liver', 'Muscle', 'Skin', 'Testis']
    if cutoff_type == 'per':
        top_coln = '%.d_percentile_counts'%(100*cutoff)
        bot_coln = '%.d_percentile_counts'%(100-100*cutoff)
        top_frac_coln = '%.d_percentile_frac'%(100*cutoff)
        bot_frac_coln = '%.d_percentile_frac'%(100-100*cutoff)
    elif cutoff_type == 'value':
        top_coln = '%s_greater_%s_counts'%(variable, cutoff[0])
        bot_coln = '%s_smaller_%s_counts'%(variable, cutoff[1])
        top_frac_coln = '%s_greater_%s_frac'%(variable, cutoff[0])
        bot_frac_coln = '%s_smaller_%s_frac'%(variable, cutoff[1])     
    for tissue in tissues:
        df_tissue = df_CM[(df_CM['Tissue'] == tissue) & (df_CM['Sample.Type'] == st) & (np.isfinite(df_CM[variable]))]
        if df_tissue.empty == False:
            if cutoff_type == 'per':
                tissue_top = df_tissue[df_tissue[variable] >= df_tissue[variable].quantile(cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                tissue_bottom = df_tissue[df_tissue[variable] <= df_tissue[variable].quantile(1-cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            elif cutoff_type == 'value':
                if p_coln == None:
                    tissue_top = df_tissue[df_tissue[variable] >= cutoff[0]].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                    tissue_bottom = df_tissue[df_tissue[variable] <= cutoff[1]].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                else:
                    tissue_top = df_tissue[(df_tissue[variable]>=cutoff[0]) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                    tissue_bottom = df_tissue[(df_tissue[variable]<= cutoff[1]) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            df_top = pd.DataFrame(tissue_top.reset_index()).rename(columns = {'N. furzeri Final Symbol':top_coln})
            df_bottom = pd.DataFrame(tissue_bottom.reset_index()).rename(columns = {'N. furzeri Final Symbol':bot_coln})
            tissue_all = df_tissue.groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            df_all = pd.DataFrame(tissue_all.reset_index()).rename(columns = {'N. furzeri Final Symbol':'all counts'})
            df_summary = pd.merge(df_top, df_bottom, how = 'outer', left_on = ['name', 'Tissue'], right_on = ['name', 'Tissue'])
            df_summary = pd.merge(df_summary, df_all, how = 'outer', left_on = ['name', 'Tissue'], right_on = ['name', 'Tissue'])
            df_summary[top_frac_coln] = df_summary[top_coln]*1.0/df_summary['all counts'].sum()
            df_summary[bot_frac_coln] = df_summary[bot_coln]*1.0/df_summary['all counts'].sum()
            frames.append(df_summary)
    df_out = pd.concat(frames, ignore_index = True)
    df_out.fillna(value = 0, inplace = True)
    df_out['rest'] = df_out['all counts'] - df_out[top_coln] - df_out[bot_coln]
    df_out['rest_frac'] = 1 - df_summary[top_frac_coln] - df_summary[bot_frac_coln]
    return df_out, top_coln, bot_coln


# do accounting on fraction of proteins (with respect to total detected in that tissue) that are in top or bottom 10 percentile
# or meet certain cutoff in each tissue, then figure out where they fall
def tissue_CM_frac(df_CM, st, variable, cutoff, sign = 'greater', cutoff_type='per', p_coln = None, p_cutoff =0.05):
    tissues = df_CM['Tissue'].unique()
    frames = []
    if cutoff_type  == 'per':
        out_coln = '>=%.d_percentile_counts'%(100*cutoff)
    elif cutoff_type == 'value':
        out_coln = '%s_>=%s_counts'%(variable, cutoff)
    elif cutoff_type == 'all':
        out_coln = 'all'
    # note the fraction here is calculated as the frac of variable that appear in a particular
    out_frac_coln = out_coln.replace('counts','frac')
    for tissue in tissues:
        df_tissue = df_CM[(df_CM['Tissue'] == tissue) & (df_CM['Sample.Type'] == st) & (np.isfinite(df_CM[variable]))]
        if df_tissue.empty == False:
            if sign == 'greater':
                if cutoff_type == 'per':
                    tissue_top = df_tissue[df_tissue[variable]>=df_tissue[variable].quantile(cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                elif cutoff_type == 'value':
                    if p_coln == None:
                        tissue_top = df_tissue[df_tissue[variable]>=cutoff].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                    else:
                        tissue_top = df_tissue[(df_tissue[variable]>=cutoff) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            elif sign == 'smaller':
                if cutoff_type == 'per':
                    tissue_top = df_tissue[df_tissue[variable]<= df_tissue[variable].quantile(cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
                elif cutoff_type == 'value':
                    if p_coln == None:
                        tissue_top = df_tissue[df_tissue[variable]<=cutoff].groupby(['name', 'Tissue'])['N. furzeri symbol'].count()
                    else:
                        tissue_top = df_tissue[(df_tissue[variable]<=cutoff) & (df_tissue[p_coln]<=p_cutoff)].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            else: # this is case to just do accounting on the entire table with no filtering what so ever
                tissue_top = df_tissue.groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            df_top = pd.DataFrame(tissue_top.reset_index()).rename(columns = {'N. furzeri Final Symbol':out_coln})
            df_top[out_frac_coln] = df_top[out_coln]*1.0/df_top[out_coln].sum()  
            if p_coln == None:
                tissue_all = df_tissue.groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            else:
                tissue_all = df_tissue[df_tissue[p_coln]<=p_cutoff].groupby(['name', 'Tissue'])['N. furzeri Final Symbol'].count()
            df_all = pd.DataFrame(tissue_all.reset_index()).rename(columns = {'N. furzeri Final Symbol':'all counts'})
            df_summary = pd.merge(df_top, df_all, how = 'outer', left_on = ['name', 'Tissue'], right_on = ['name', 'Tissue'])
            df_summary[out_frac_coln.replace('frac','frac_ofdetected')] = df_summary[out_coln]*1.0/df_summary['all counts'].sum()            
            frames.append(df_summary)
    df_out = pd.concat(frames, ignore_index = True)
    df_out.fillna(value = 0, inplace = True)
    return df_out, out_coln, out_frac_coln


# plot a modified windrose plot to visualize the cellular localization for a compartment across different tissues
def modified_windrose(df_summary, tissues, colns, compartment, ax = 0, cm_tissues = None, cm_rest = None, cm_bot = None, cm_top = None):
    # here colns has to specify  the name of [top_col, bot_col, rest_col]
    if cm_tissues == None: cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
    if cm_rest == None: cm_rest = 'gray'
    if cm_bot == None: cm_bot = 'lightskyblue'
    if cm_top == None: cm_top = 'gold'
    df_cmi = df_summary[df_summary['name'] == compartment]
    tissues_unique = df_cmi['Tissue'].unique()
    tissues = [i for i in tissues if i in tissues_unique]
    if df_cmi.empty == True:
        return None
    else:
        top_col, bot_col,rest_col = colns
        if ax == 0: 
            ax = plt.subplot(projection = 'polar')
        N = len(tissues)
        theta = np.arange(0.0, 2*np.pi, 2*np.pi/N)
        # width represent the width of each slicE which should be evenly distribute in 2pi/N
        bot_tissue = np.asarray([10]*N); bot_rest = np.asarray([20]*N); # bot_top and bot_bot will be dynamic depending on the # of proteins
        width = [2.*np.pi/N]*N
        # the radii should represent the # of proteins in each category
        radii_bot = np.zeros(N)
        radii_top = np.zeros(N)
        radii_rest = np.zeros(N)
        for i, tissue in enumerate(tissues):
            bot_count = df_cmi[bot_col][df_cmi['Tissue'] == tissue].tolist()
            top_count = df_cmi[top_col][df_cmi['Tissue'] == tissue].tolist()
            rest_count = df_cmi[rest_col][df_cmi['Tissue'] == tissue].tolist()
            radii_bot[i] = math.sqrt(bot_count[0]) if bot_count != [] else 0
            radii_top[i] = math.sqrt(top_count[0]) if top_count != [] else 0
            radii_rest[i] = math.sqrt(rest_count[0]) if rest_count != [] else 0
        # plot all the bar charts in polar coordinates
        bars_tissue = ax.bar(theta, bot_tissue, width = width, bottom = bot_tissue, \
                             edgecolor = 'black')
        bars_rest = ax.bar(theta, radii_rest, width=width, bottom = bot_rest, \
                           edgecolor = 'black', color = cm_rest, alpha = 1)
        bars_bot = ax.bar(theta, radii_bot, width=width, bottom = radii_rest + bot_rest, \
                          edgecolor = 'black', color = cm_bot, alpha = 1)
        bars_top = ax.bar(theta, radii_top, width=width, bottom = radii_rest + radii_bot + bot_rest, \
                          edgecolor = 'black', color = cm_top, alpha = 1)
        for i, r, bar_tissue in zip (range(N), bot_tissue, bars_tissue):
            bar_tissue.set_facecolor(cm_tissues[i])
            bar_tissue.set_alpha(1)
    #    # label each tissue
    #    for i, phi in enumerate(theta):
    #        ax.annotate(tissues[i], xy = (phi, 100))
        plt.axis('off')
        return ax
    
    
# plot a modified windrose plot to visualize the cellular localization for a compartment across different tissues
def windrose_varyslice(df_summary, tissues, coln, compartment, ax = 0, cm_tissues = None):
    # here coln is one variable only!!!
    # useful references: https://matplotlib.org/3.2.1/gallery/pie_and_polar_charts/nested_pie.html#sphx-glr-gallery-pie-and-polar-charts-nested-pie-py
    if cm_tissues == None: cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
    df_cmi = df_summary[df_summary['name'] == compartment]
    tissues_unique = df_cmi['Tissue'].unique()
    tissues = [i for i in tissues if i in tissues_unique]
    if df_cmi.empty == True:
        return None
    else:
        if ax == 0: 
            ax = plt.subplot(projection = 'polar')
        N = len(tissues)
        theta = np.linspace(0, 2*np.pi,N)#np.zeros(N) #the angle here is dynamic depending on the # of proteins
        # width represent the width of each slicE which should be determined by  # of proteins
        width = np.zeros(N) #np.zeros(N)#[2.*np.pi/N]*N
        tot = df_cmi[coln].sum()
        for i, tissue in enumerate(tissues):
            width[i] = 2*np.pi * df_cmi[coln][df_cmi['Tissue'] == tissue].values[0]/tot
            theta[i] = width[:i].sum()
        # plot all the bar charts in polar coordinates
        size = 0.5
        bars_tissue = ax.bar(x = theta, width = width, bottom = 1-size, height = size, edgecolor = 'black', align='edge')
#        if tot <10:
#            ax.annotate('%s'%int(tot), xy=(0, 0), xytext = (0.45, 0.45), fontsize=48, textcoords='axes fraction')
#        elif tot <100:
#            ax.annotate('%s'%int(tot), xy=(0, 0), xytext = (0.4, 0.45), fontsize=48, textcoords='axes fraction')
#        else:
#            ax.annotate('%s'%int(tot), xy=(0, 0), xytext = (0.35, 0.45), fontsize=48, textcoords='axes fraction')
        avg_frac = 100*df_cmi[coln].mean()
        if avg_frac>10:
            ax.annotate('%.1f%%'%avg_frac, xy=(0, 0), xytext = (0.3, 0.45), fontsize=32, textcoords='axes fraction')
        else:
            ax.annotate('%.1f%%'%avg_frac, xy=(0, 0), xytext = (0.35, 0.45), fontsize=32, textcoords='axes fraction')
        for i, bar_tissue in zip (range(N), bars_tissue):
            bar_tissue.set_facecolor(cm_tissues[i])
            bar_tissue.set_alpha(1)
        plt.axis('off')
        return ax


def cm_windrose(df_summary, tissues, colns, compartments, fname, fpath, windrose_type ='norm', cm_tissues = None, cm_rest = None, cm_bot = None, cm_top = None):    
    # use sharex = True and sharey = True to set the same absolute size
#    fig, axes = plt.subplots(nrows = 4, ncols = 5, figsize = (24, 24), subplot_kw={'projection':'polar'})
    fig, axes = plt.subplots(nrows = 4, ncols = 5, figsize = (24, 24), subplot_kw={'projection':'polar'}, sharex = True, sharey = True)
    ax_list = axes.flat
    for i, cm in enumerate(compartments): 
        if cm not in df_summary['name'].values:
            continue
        if windrose_type == 'norm':
            ax_i = modified_windrose(df_summary, tissues, colns, cm, ax = ax_list[i], \
                                 cm_tissues = cm_tissues, cm_rest = cm_rest, cm_bot = cm_bot, cm_top = cm_top)
        elif windrose_type == 'varyslice':
            ax_i = windrose_varyslice(df_summary, tissues, colns, cm, ax = ax_list[i], cm_tissues = cm_tissues)
        if ax_i != None:
            ax_i.set_title(cm, fontsize=32)
            ax_i.axis('off')
    fig.suptitle(fname.replace('.pdf',''))
    fig.savefig(os.path.join(fpath, fname))
    plt.close(fig)
    return None


def modified_pie(slice_sizes, radius, slice_colors = None, ax = None):
    if ax == None:
        fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'polar'})
    N = len(slice_sizes)
    if slice_colors == None:
        cmap = mpl.cm.get_cmap('viridis', N)
        slice_colors = [cmap(i) for i in np.asarray(slice_sizes)/sum(slice_sizes)]
    # width should be the same since we are generating a pie chart like circle
    width = np.asarray(slice_sizes)*1./sum(slice_sizes)*2.*np.pi
    # distribute the angles based on the slice sizes across 2pi
    theta_sum = 0
    theta = []
    for i in width:
        theta.append(i/2. + theta_sum)
        theta_sum += i
    theta = np.asarray(theta)
    # the radii should represent the # of proteins in each category
    radii = np.ones(N)*radius
    # plot all the bar charts in polar coordinates
    slices = ax.bar(theta, radii, width = width, bottom = 0, edgecolor = 'black')
    for i, r, bar in zip (range(N), radii, slices):
        bar.set_facecolor(slice_colors[i])
        bar.set_alpha(1)
#    ax.axis('off')
    return ax


# plot a modified pizza slice-like circle to visualize the cellular localization for Increased/Decrease AGG/Prop for a compartment across different tissues
def pizza_spike(df_summary, tissues, col, compartment, ax = None, cm_tissues = None, count = False):
    if cm_tissues == None: cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
    df_cmi = df_summary[df_summary['name'] == compartment]
    tissues_unique = df_cmi['Tissue'].unique()
    tissues = [i for i in tissues if i in tissues_unique]
    if df_cmi.empty == True:
        return None
    else:
#        if sign == 'greater': 
#            top_col = '>=%.d percentile counts'%(100*per) if count == True else '>=%.d_frac'%(100*per)
#        elif sign == 'smaller':
#            top_col = '<=%.d percentile counts'%(100*per) if count == True else '<=%.d_frac'%(100*per)
        if ax == None: 
            ax = plt.subplot(projection = 'polar')
        N = len(tissues)
        theta = np.arange(0.0, 2*np.pi, 2*np.pi/N)
        # width represent the width of each slicE which should be evenly distribute in 2pi/N
        bot_rest = np.asarray([0.5]*N) if count == True else np.asarray([0.001]*N) # bot_top and bot_bot will be dynamic depending on the # of proteins
        width = [2.*np.pi/N]*N
        # the radii should represent the # of proteins in each category
        radii_top = np.zeros(N)
        for i, tissue in enumerate(tissues):
            top_count = df_cmi[col][df_cmi['Tissue'] == tissue].tolist()
            radii_top[i] = math.sqrt(top_count[0]) if top_count != [] else 0
        # plot all the bar charts in polar coordinates
        bars_top = ax.bar(theta, radii_top, width=width, bottom = bot_rest, edgecolor = 'black', alpha = 1)
        for i, bar_tissue in zip (range(N), bars_top):
            bar_tissue.set_facecolor(cm_tissues[i])
            bar_tissue.set_alpha(1)
    #    # label each tissue
    #    for i, phi in enumerate(theta):
    #        ax.annotate(tissues[i], xy = (phi, 100))
        plt.axis('off')
        return ax


def cm_pizza_spike(df_summary, tissues, col, compartments, fname, fpath, cm_tissues = None, count = False):    
    # use sharex = True and sharey = True to set the same absolute size
#    fig, axes = plt.subplots(nrows = 4, ncols = 5, figsize = (24, 24), subplot_kw={'projection':'polar'})
    fig, axes = plt.subplots(nrows = 4, ncols = 5, figsize = (24, 24), subplot_kw={'projection':'polar'}, sharex = True, sharey = True)
    ax_list = axes.flat
    for i, cm in enumerate(compartments): 
        ax_i = pizza_spike(df_summary, tissues, col, cm, ax = ax_list[i], cm_tissues = cm_tissues, count = count)
        if ax_i != None:
            ax_i.set_title(cm)
            ax_i.axis('off')
    fig.savefig(os.path.join(fpath, fname)) if count == True else fig.savefig(os.path.join(fpath, fname.replace('.pdf','_frac.pdf')))
    return None

def new_rank(df, fc, pval):
    df['rankstats'] = df[fc]* -np.log10(df[pval]) 
    return df


def cm_grid(df_summary, tissues, compartments, per, variable, fpath, fname, cm_tissues = None):
    if cm_tissues == None: cm_tissues = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
    fig, axes = plt.subplots(nrows = len(compartments), ncols = len(tissues), figsize = (7,20), \
                             subplot_kw={'projection':'polar'}, sharex = 'row', sharey = 'col')
    # the plot will resemble a scatter plot where x is tissue and y is compartments or vice versa but instead plotted through subplots
    ax_list = axes.flat
    i = 0
    top_col = '%.d percentile counts'%(100*per)
    bot_col = '%.d percentile counts'%(100-100*per)
    rest_col = 'rest'
    for cm_i in compartments:
        df_cmi = df_summary[df_summary['name'] == cm_i]
        for tissue in tissues:
            ax_i = ax_list[i]
            top = df_cmi[top_col][df_cmi['Tissue'] == tissue].tolist()[0]
            bot = df_cmi[bot_col][df_cmi['Tissue'] == tissue].tolist()[0]
            rest = df_cmi[rest_col][df_cmi['Tissue'] == tissue].tolist()[0]
            modified_pie([math.sqrt(top), math.sqrt(bot), math.sqrt(rest)], math.sqrt(top + bot+ rest), ax = ax_i, slice_colors = ['gold', 'lightskyblue', 'gray'])
            ax_i_boolean = False
            if i % len(tissues) == 0:
                ax_i.set_ylabel(cm_i)
                # hide the x and y axis ticks but not the labels
                ax_i.get_xaxis().set_ticks([])
                ax_i.get_yaxis().set_ticks([])
                ax_i_boolean = True
            if cm_i == compartments[-1]:
                ax_i.set_xlabel(tissue)
                # hide the x and y axis ticks but not the labels
                ax_i.get_xaxis().set_ticks([])
                ax_i.get_yaxis().set_ticks([])
                ax_i_boolean = True
            if ax_i_boolean == False:
                ax_i.axis('off')
#            ax_i.get_yaxis().set_visible(False)
#            ax_i.axis('off')
            i += 1
#    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#    plt.xlabel('brain liver gut')
#    plt.tight_layout()
    fig.savefig(os.path.join(fpath, fname))
    return None

###### SECTION 1: OBTAIN CELLULAR LOCALIAZTION ON TABLE OF INTEREST (TOI) ####
## updated on 2020.05.04    
fpath = 'CellularLocalization'
datapath = ''
CM_database = 'killifish_human_CM.csv'
TOI = 'kf_combined_prop.csv'
df_kf_localization = pd.read_csv(os.path.join ('CellularLocalization/Uniprot', CM_database)) 
df_MS = pd.read_csv(os.path.join(datapath, TOI))
df_MS.drop(df_MS.columns[76:126],axis = 1, inplace=True)
# merge the killifish localization data with the TMT-MS data
df_CM = df_MS.merge(df_kf_localization, how = 'left', left_on = ['N. furzeri Protein Id', 'N. furzeri Gene Id'], \
                    right_on = ['Kf Id2', 'Kf Id'])
df_CM.to_csv(os.path.join(fpath, TOI.replace('.csv', '_CM.csv')),index=False)
############################## END OF SECTION 1 ##############################

######## SECTION 3: OBTAIN CELLULAR LOCALIAZTION STATS BASED ON PERCENTILE ######    
## PLOT THEM IN A PIE CHART FORMAT FOR EACH TISSUE AND EACH COMPARTMENT #
# updated on 2020.05.04
cutoff = None #.75
p_cutoff = .05#.05
age1 = 'Young'; age2 = 'Old'; age_comp ='%sv%s'%(age2[0], age1[0])
st = 'AGG'
fout_name = '%s_%s_Windrose_scale.pdf'%(st, age_comp) if cutoff == None else '%sZ%s_%s_Windrose_scale.pdf'%(st, cutoff, age_comp)
if st == 'Prop':
    fc = '%s_prop_logFC'%age_comp
    pval = '%s_prop_pval'%age_comp
else:
    fc = '%s_logFC'%age_comp
    pval = '%s_pval'%age_comp
# initiaze variables    
tissues = ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']
tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
compartments = ['Golgi', 'cytoplasm', 'cytosol', 'endoplasmic reticulum', 'extracellular region', \
               'lysosome', 'mitochondrial ribosome', 'mitochondrion', 'nuclear envelope', 'nuclear pore',\
               'nucleolus', 'nuclear speck', 'nucleoplasm', 'nucleus', 'plasma membrane', \
               'peroxisome', 'proteasome', 'ribosome', 'secretory granule', 'splicesome']
## load Cellular Compartment data
fpath = 'CellularLocalization'
fname_CM = 'kf_combined_prop_CM.csv'#'test_CM.csv'
fpath_out = 'CellularLocalization/ByFCzscore'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
df_CM = pd.read_csv(os.path.join(fpath, fname_CM))
######## part 1. find out what the z-score>0.75 aggregating protein and their cellular localization for each tissue
v_cutoff = cutoff if cutoff != None else 0
for sign in ['smaller', 'greater']:
    direction = 'Pos' if sign == 'greater' else 'Neg'
    if st == 'Prop':
        df_CM_summary, top_col, topfrac_col = tissue_CM_frac(df_CM, 'AGG', fc, v_cutoff, sign = sign, cutoff_type ='value', p_coln = pval, p_cutoff = p_cutoff)
    else:
        df_CM_summary, top_col, topfrac_col = tissue_CM_frac(df_CM, st, fc, v_cutoff, sign = sign, cutoff_type ='value', p_coln = pval, p_cutoff = p_cutoff)
    df_CM_out = df_CM_summary.pivot(index='name', columns='Tissue', values=list(set(df_CM_summary.columns) - set(['name','Tissue'])))
    df_CM_out.to_csv(os.path.join(fpath_out, fout_name.replace('Windrose_scale.pdf','%sFrac.csv'%direction)))
    cm_windrose(df_CM_summary, tissues, top_col.replace('counts','frac'), compartments, fout_name.replace('Windrose_scale.pdf','%sFrac_Donut.pdf'%direction), fpath_out, windrose_type ='varyslice')
################################ END OF SECTION 3 ##############################    
    

######## SECTION 4: OBTAIN CELLULAR LOCALIAZTION STATS FOR DETECTED PROTEINS ######
## PLOT THEM IN A WINDROSE STYLE PLOT TO SHOW DIFFERENCES AMONG DIFFERENT TISSUES (PER COMPARTMENT)
# updated on 2020.05.04
age = 'Young'; st = 'TL'; cutoff = None
if cutoff == None: fout_name = '%s_All_Windrose_scale.pdf'%(st)
var = 'log2_%s_prop'%age if st == 'Prop' else 'log2_%s'%age
# initiaze variables    
tissues = ['Brain','Gut','Heart','Liver','Muscle','Skin','Testis']
tissues_colors = ['#d53e4f', '#fc8d59', '#fee08b', '#ffffbf', '#e6f598', '#99d594', '#3288bd']
compartments = ['Golgi', 'cytoplasm', 'cytosol', 'endoplasmic reticulum', 'extracellular region', \
              'lysosome', 'mitochondrial ribosome', 'mitochondrion', 'nuclear envelope', 'nuclear pore',\
              'nucleolus', 'nuclear speck', 'nucleoplasm', 'nucleus', 'plasma membrane', \
              'peroxisome', 'proteasome', 'ribosome', 'secretory granule', 'splicesome']
## load Cellular Compartment data
fpath = 'CellularLocalization'
fname_CM = 'kf_combined_prop_CM.csv'
fpath_out = 'CellularLocalization/ByExp'
if not os.path.exists(fpath_out): os.makedirs(fpath_out)
df_CM = pd.read_csv(os.path.join(fpath, fname_CM))
######## part 1. find out the cellular localization for detected proteome in each tissue
v_cutoff = cutoff if cutoff != None else 0
if st == 'Prop':
    df_CM_summary, top_col, topfrac_col = tissue_CM_frac(df_CM, 'AGG', var, v_cutoff, sign = 'all')
else:
    df_CM_summary, top_col, topfrac_col = tissue_CM_frac(df_CM, st, var, v_cutoff, sign = 'all')
df_CM_out = df_CM_summary.pivot(index='name', columns='Tissue', values=list(set(df_CM_summary.columns) - set(['name','Tissue'])))
df_CM_out.to_csv(os.path.join(fpath_out, fout_name.replace('Windrose_scale.pdf','Frac.csv')))
cm_windrose(df_CM_summary, tissues, top_col.replace('counts','frac'), compartments, fout_name.replace('Windrose_scale.pdf','Frac_Donut.pdf'), '%s/ByExp'%fpath, windrose_type ='varyslice')