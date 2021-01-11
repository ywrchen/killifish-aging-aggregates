#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 16:32:05 2020

@author: yiwenchen
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import pandas as pd
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 1)
sns.set(style="white")


def dot_plot(df, x, y, hue=None, ax=None, order=None):
    sns.stripplot(x=x, y=y, hue=hue, data=df, order=order,\
                  jitter=True, dodge=True, ax=None,\
                  edgecolor='k', linewidth=1, size=10)
    labels = [e.get_text() for e in plt.gca().get_xticklabels()]
    ticks = plt.gca().get_xticks()
    w = 0.12
    for i, strain in enumerate(labels):
        idx = labels.index(strain)
        for j, cond in enumerate(df[hue].unique()):
            adj = 2*w*(j-len(df[hue].unique())+2)
            plt.hlines(df[(df[x] == strain)&(df[hue] == cond)][y].mean(), ticks[idx]-w+adj, ticks[idx]+w+adj, color='k')
    ax.set(ylim=(0,1))
    return ax

    
fname_data = '20200821_YC_TwoColor_CPQuant.csv'
fpath_data = 'BrainFig'

df_data = pd.read_csv(os.path.join(fpath_data, fname_data))
convert_dict = {'%wGFPpucta.per.GFPpositive': float, 
                '%wGFPpuncta.per.cell': float} 
df_data = df_data.astype(convert_dict) 

#####  Section 1: main figure 
fig1, ax1 = plt.subplots(1,1, figsize=(8,4))
strain_order = ['YDJ5500', 'YDJ5501', 'YDJ6615', 'YDJ5530', 'YDJ5531', 'YDJ6616']
df_fig1 = df_data[df_data['Strain'].isin(strain_order)]
ax = dot_plot(df_fig1, x='Strain', y='%wGFPpucta.per.GFPpositive', hue = 'Condition', \
              order=strain_order, ax=ax1)   
fig1.savefig(os.path.join(fpath_data,'TwoColor_main1.pdf'))

##### Section 2: sup figure (polyQ)
fig2, ax2 = plt.subplots(1,1, figsize=(4,4))
strain_order = ['YDJ6617', 'YDJ6618', 'YDJ6619']
df_fig2 = df_data[df_data['Strain'].isin(strain_order)]
ax2 = dot_plot(df_fig2, x='Strain', y='%wGFPpucta.per.GFPpositive', hue = 'Condition', \
              order=strain_order, ax=ax2)   
fig2.savefig(os.path.join(fpath_data,'TwoColor_main2.pdf'))

##### Section 2: sup figure (polyQ)
fig3, ax3 = plt.subplots(1,1, figsize=(2,4))
strain_order = ['YDJ5502']
df_fig3 = df_data[df_data['Strain'].isin(strain_order)]
ax3 = dot_plot(df_fig3, x='Strain', y='%wGFPpucta.per.GFPpositive', hue = 'Condition', \
              order=strain_order, ax=ax3)   
fig3.savefig(os.path.join(fpath_data,'TwoColor_main3.pdf'))
