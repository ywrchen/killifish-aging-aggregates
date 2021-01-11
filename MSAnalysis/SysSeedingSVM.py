#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# for details read https://scikit-learn.org/stable/auto_examples/svm/plot_separating_hyperplane.html
# get the hyperplane https://scikit-learn.org/0.17/auto_examples/svm/plot_separating_hyperplane.html
# https://medium.com/deep-math-machine-learning-ai/chapter-3-support-vector-machine-with-math-47d6193c82be

"""
Created on Thu Dec  3 09:17:14 2020

@author: yiwenchen
"""

import pandas as pd
import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import itertools
import random
from scipy import stats
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV

mpl.rcParams['pdf.fonttype'] = 42
sns.set(font_scale = 4)
sns.set_style("white")

class MidpointNormalize(Normalize):

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# # ################### Section 1: Build 2 parameter model and define hyperplane (support vector machine) ##################
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# features = ['OvY_logFC', 'OvY_prop_logFC','log2_Young_TL',\
#             'AA Length', 'FImaxrun','Longest Disorder Stretch',\
#             'FracExpanding','FracDisoPromoting', 'Disordered Frac','Frac_FInumaa',\
#             'PAPAprop','NLLR','MWScore','Frac_QN',\
#             'pI', 'FracNeg','Frac_pos','FracCharged',\
#             'delta', 'kappa','Omega','MeanNetCharge','uversky_hydropathy','Frac_neu',\
#             'Frac_aromatic','Frac_aliphatic']
#             # 'mean Ka/Ks','mean Ka','mean ks'] coverage for these metrices are poor so excluded in the analysis
# class_coln = 'Aggregate'
# X = df_data[features].values
# Y = df_data[class_coln]
# ## visualize distribution of features
# #df_data[features].plot(kind='box', subplots=True, layout=(5,5), sharex=False, sharey=False, figsize=(30,30))
# #df_data[features].hist(bins=10, figsize=(30,30))
# ## visaulize features that might be correlated
# #scatter = pd.plotting.scatter_matrix (df_data[features], c=Y.values, marker = 'o', s=40, figsize=(20,20))
# X = StandardScaler().fit_transform(X)
# df_scaled = pd.DataFrame(X, index=df_data.index, columns = features)
# # iterate through pair of 2 parameters and identify ones that best separate two classes with linear kernel
# # feature_pairs = list(itertools.combinations(features,r=2))
# feature_pairs = list(itertools.combinations_with_replacement(features, r=2)) 
# # use support vector classifier to find the hyperplane that maximize margin
# clf = svm.SVC(kernel='linear', C=1)
# #X = MinMaxScaler().fit_transform(X)
# #trainedsvm = svm_model.fit(X_train,Y_train)
# clf_out = []
# n_test = 50
# random.seed(5)
# r = random.sample(range(1000),n_test)
# for f_pair in feature_pairs:
#     var_x, var_y = f_pair
#     df_pair = df_scaled[[var_x,var_y]].dropna(axis=0,how='any')
#     X_pair = df_pair.values
#     Y_pair = df_data.loc[df_pair.index,class_coln]
#     clf.fit(X_pair, Y_pair)
#     Y_pred = clf.predict(X_pair)
#     accuracy_all = metrics.accuracy_score(Y_pair, Y_pred)
#     clf_i = [accuracy_all]
#     for i in range(n_test):
#         X_train, X_test, Y_train, Y_test = train_test_split(X_pair, Y_pair, test_size = 0.2,random_state=r[i])
#         clf.fit(X_train, Y_train)
#         Y_test_pred = clf.predict(X_test)
#         accuracy = metrics.accuracy_score(Y_test, Y_test_pred)
#         clf_i.append(accuracy) 
#     clf_out.append([var_x,var_y]+clf_i)
# df_clf = pd.DataFrame(clf_out, columns = ['var_x','var_y','Score_All']+['Score%s'%(i+1) for i in range(n_test)])
# df_clf['Score_100mean'] = df_clf[['Score%s'%(i+1) for i in range(n_test)]].mean(axis=1)*100
# df_clf.to_csv(os.path.join(out_dir, 'SysSeedingLinearSVC.csv'), index=False)      
# # ################### End of Section 1 ###################  

  
# ###################  Section 2: Visualize some of the best SVM classifier (scaled) ################### 
# sns.set(font_scale = 4)
# var_x = 'OvY_prop_logFC'; var_y = 'delta' ;C=1; gamma=1
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# class_coln = 'Aggregate'
# X_raw = df_data[[var_x,var_y]].dropna(axis=0,how='any')
# Y = df_data.loc[X_raw.index, class_coln]
# frac_mapped = X_raw.shape[0]*1./df_data.shape[0]
# # standarize the input data
# X = StandardScaler().fit_transform(X_raw.values)
# df_scaled = pd.DataFrame(X, index=X_raw.index, columns = [var_x,var_y])
# df_scaled['Tissue'] = df_data.loc[X_raw.index, 'Tissue']
# df_scaled[class_coln] = df_data.loc[X_raw.index, class_coln]
# # use support vector classifier to find the hyperplane that maximize margin
# clf = svm.SVC(kernel='linear', C=C, gamma=gamma)
# clf.fit(X, Y)
# Y_pred = clf.predict(X)
# accuracy_all = metrics.accuracy_score(Y, Y_pred)
# # plot original data
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na', 
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_scaled , ax = ax, \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# # plot the decision function
# xlim = ax.get_xlim()
# ylim = ax.get_ylim()
# # create grid to evaluate model
# xx = np.linspace(xlim[0], xlim[1], 30)
# yy = np.linspace(ylim[0], ylim[1], 30)
# YY, XX = np.meshgrid(yy, xx)
# xy = np.vstack([XX.ravel(), YY.ravel()]).T
# Z = clf.decision_function(xy).reshape(XX.shape)
# # plot decision boundary and margins
# ax.contour(XX, YY, Z, colors='k', levels=[0], alpha=0.5, linestyle=['-'])
# ax.text(0.05,0.9, 'accuracy = %.2f%%\ncoverage = %.2f%%'%(accuracy_all*100,frac_mapped*100), transform=ax.transAxes)
# #ax.contour(XX, YY, Z, colors='k', levels=[-1,0,1], alpha=0.5, linestyle=['--','-','--'])
# #    # plot the support vector
# #    sns.scatterplot(clf.support_vectors_[:,0], clf.support_vectors_[:,1], linewidth=2, facecolors='none', s=5000)
# # add legend etc.
# ax.legend(markerscale=6, bbox_to_anchor=(1.4, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(out_dir, 'SysSeeding_%sv%s_scaled.pdf'%(var_y, var_x)))
# print('Accuracy of SVM classifier with %s and %s : %.2f'%(var_x, var_y,accuracy_all))
# plt.close()
# #################### End of Section 2 ################### 


# ##################  Section 3: Visualize some of the best SVM classifer (un-transformed coordinate) ################### 
# var_x = 'OvY_prop_logFC'; var_y = 'delta';C=1; gamma=1 #var_x = 'OvY_prop_logFC', 'log2_Young_TL'
# sns.set(font_scale = 4)
# sns.set_style("white")
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# class_coln = 'Aggregate'
# X_raw = df_data[[var_x,var_y]].dropna(axis=0,how='any')
# Y = df_data.loc[X_raw.index, class_coln]
# frac_mapped = X_raw.shape[0]*1./df_data.shape[0]
# ### standarize the input data
# scal = StandardScaler()
# X = scal.fit_transform(X_raw.values)
# df_scaled = pd.DataFrame(X, index=X_raw.index, columns = [var_x,var_y])
# df_scaled['Tissue'] = df_data.loc[X_raw.index, 'Tissue']
# df_scaled[class_coln] = df_data.loc[X_raw.index, class_coln]
# # use support vector classifier to find the hyperplane that maximize margin
# clf = svm.SVC(kernel='linear', C=C, gamma=gamma)
# clf.fit(X, Y)
# Y_pred = clf.predict(X)
# accuracy_all = metrics.accuracy_score(Y, Y_pred)
# # plot original data
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na', 
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_data , ax = ax, \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# # get the separating hyperplane in original coordinate
# w = clf.coef_[0]
# a = -w[0]/w[1]
# xx = X[:,0]
# yy = a*xx - (clf.intercept_[0])/w[1]
# raw_plane = scal.inverse_transform([(i,j) for i, j in zip(xx,yy)])
# # plot the separating hyperplane/lines
# ax.plot(raw_plane[:,0], raw_plane[:,1], 'k-')
# ax.text(0.05,0.9, 'accuracy = %.2f%%\ncoverage = %.2f%%'%(accuracy_all*100,frac_mapped*100), transform=ax.transAxes)
# # add legend etc.
# ax.legend(markerscale=6, bbox_to_anchor=(1.4, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(out_dir, 'SysSeeding_%sv%s.pdf'%(var_y, var_x)))
# print('Accuracy of SVM classifier with %s and %s : %.2f'%(var_x, var_y,accuracy_all))
# plt.close()
#################### End of Section 3 ###################

# ##################  Section 4: Visualize some of the best SVM classifer for output (scaled, explicit hyperplane) ###################
# var_x = 'OvY_prop_logFC';var_y = 'delta' #var_x = 'OvY_prop_logFC'
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# class_coln = 'Aggregate'
# X_raw = df_data[[var_x,var_y]].dropna(axis=0,how='any')
# Y = df_data.loc[X_raw.index, class_coln]
# frac_mapped = X_raw.shape[0]*1./df_data.shape[0]
# ### standarize the input data
# scal = StandardScaler()
# X = scal.fit_transform(X_raw.values)
# df_scaled = pd.DataFrame(X, index=X_raw.index, columns = [var_x,var_y])
# df_scaled['Tissue'] = df_data.loc[X_raw.index, 'Tissue']
# df_scaled[class_coln] = df_data.loc[X_raw.index, class_coln]
# # use support vector classifier to find the hyperplane that maximize margin
# clf = svm.SVC(kernel='linear', C=1)
# clf.fit(X, Y)
# Y_pred = clf.predict(X)
# accuracy_all = metrics.accuracy_score(Y, Y_pred)
# # plot original data
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na',
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#               style = 'Tissue', data = df_scaled, ax = ax, \
#               palette = aggregate_colors, hue_order = aggregate_order, s=5000)
# # plot the decision function
# xlim = ax.get_xlim()
# ylim = ax.get_ylim()
# # create grid to evaluate model
# xx_plane = np.linspace(xlim[0], xlim[1], 30)
# yy_plane = np.linspace(ylim[0], ylim[1], 30)
# scaled_plane = scal.fit_transform([(i,j) for i, j in zip(xx_plane,yy_plane)])
# # get the separating hyperplane
# w = clf.coef_[0]
# a = -w[0]/w[1]
# xx_scaled = xx_plane
# yy_scaled = a*xx_plane - (clf.intercept_[0])/w[1]
# ## plot the parallels to the separating hyperplane that pass through the support vectors
# #b = clf.support_vectors_[1]
# #yy_down = a * xx_scaled + (b[1] - a * b[0])
# #b = clf.support_vectors_[-1]
# #yy_up = a * xx_scaled + (b[1] - a * b[0])
# # plot the lines
# ax.plot(xx_scaled, yy_scaled, 'k-')
# #ax.plot(xx_scaled, yy_down, 'k--')
# #ax.plot(xx_scaled, yy_up, 'k--')
# # plot the accuracy of fit
# ax.text(0.05,0.9, 'accuracy = %.2f%%\ncoverage = %.2f%%'%(accuracy_all*100,frac_mapped*100), transform=ax.transAxes)
# # plot the support vector
# # sns.scatterplot(clf.support_vectors_[:,0], clf.support_vectors_[:,1], linewidth=2, facecolors='none', s=5000)
# # add legend etc.
# ax.legend(markerscale=6, bbox_to_anchor=(1.4, 1))
# plt.tight_layout()
# fig.savefig(os.path.join(out_dir, 'SysSeeding_%sv%s_scaled.pdf'%(var_y, var_x)))
# print('Accuracy of SVM classifier with %s and %s : %.2f'%(var_x, var_y,accuracy_all))
# plt.close()
# #################### End of Section 4 ###################


# #################### Section 5: Fine tuning a model to find the best hyperplane (support vector machine) ##################
# var_x = 'OvY_prop_logFC'; var_y = 'delta' 
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# class_coln = 'Aggregate'
# X_raw = df_data[[var_x,var_y]].dropna(axis=0,how='any')
# Y = df_data.loc[X_raw.index, class_coln]
# frac_mapped = X_raw.shape[0]*1./df_data.shape[0]
# ### standarize the input data
# scal = StandardScaler()
# X = scal.fit_transform(X_raw.values)
# df_scaled = pd.DataFrame(X, index=X_raw.index, columns = [var_x,var_y])
# df_scaled['Tissue'] = df_data.loc[X_raw.index, 'Tissue']
# df_scaled[class_coln] = df_data.loc[X_raw.index, class_coln] 
# # Train classifiers, for an initial search, a logarithmic grid with basis 10 
# C_range = np.arange(0.01, 10, 0.01)
# gamma_range = np.logspace(-9, 3, 13)
# hyperparameter_pairs = list(itertools.product(C_range,gamma_range))
# clf_out = np.zeros((len(hyperparameter_pairs),3))
# for i, pair in enumerate(hyperparameter_pairs):
#     C_i, gamma_i = pair
#     clf = svm.SVC(kernel='linear', C=C_i, gamma=gamma_i)
#     clf.fit(X, Y)
#     Y_pred = clf.predict(X)
#     accuracy_i= metrics.accuracy_score(Y, Y_pred)
#     clf_out[i,0] = C_i
#     clf_out[i,1] = gamma_i
#     clf_out[i,2] = accuracy_i
# df_clf = pd.DataFrame(data=clf_out, columns=['C','gamma','accuracy'])
# score = df_clf.pivot(index='C',columns='gamma',values='accuracy')
# fig, ax = plt.subplots(1, 1, figsize = (30,20))
# sns.heatmap(score, ax=ax)
# ax.set_title('Validation accuracy')
# fig.savefig(os.path.join(out_dir,'SysSeeding_%sv%s_tuning.pdf'%(var_y, var_x)))
# #################### End of Section 5 ###################

# #################### Section 6: plot the entire detected population
# var_x = 'log2_Young_TL'; var_y = 'delta';C=1; gamma=1 #var_x = 'OvY_prop_logFC'
# sns.set(font_scale = 1)
# sns.set_style("white")
# seeding_dir = 'SysSeeding'
# out_dir = 'SysSeeding/SVM'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# df_data = df_data[(df_data['Aggregate']!='TBD') & (df_data['Hit.Type'].isin(['AGG','Prop','Prop, AGG']))]
# df_data.loc[df_data['Aggregate']=='frac agg', 'Aggregate'] = 'aggregate'
# class_coln = 'Aggregate'
# X_raw = df_data[[var_x,var_y]].dropna(axis=0,how='any')
# Y = df_data.loc[X_raw.index, class_coln]
# frac_mapped = X_raw.shape[0]*1./df_data.shape[0]
# ### standarize the input data
# scal = StandardScaler()
# X = scal.fit_transform(X_raw.values)
# df_scaled = pd.DataFrame(X, index=X_raw.index, columns = [var_x,var_y])
# df_scaled['Tissue'] = df_data.loc[X_raw.index, 'Tissue']
# df_scaled[class_coln] = df_data.loc[X_raw.index, class_coln]
# # use support vector classifier to find the hyperplane that maximize margin
# clf = svm.SVC(kernel='linear', C=C, gamma=gamma)
# clf.fit(X, Y)
# Y_pred = clf.predict(X)
# accuracy_all = metrics.accuracy_score(Y, Y_pred)
# # plot original data
# aggregate_colors = ['#4dac26', '#d01c8b', '#f1b6da'] #'#D1CCCB'
# aggregate_order = ['diffuse', 'aggregate', 'frac agg'] #'na', 
# fig, ax = plt.subplots(1, 2, figsize = (15,7), sharex=True, sharey=True)
# sns.scatterplot(x = var_x, y = var_y, hue = 'Aggregate', \
#                 style = 'Tissue', data = df_data , ax = ax[0], \
#                 palette = aggregate_colors, hue_order = aggregate_order, s=500)
# # for line in df_data.index:
# #     # x_offset = 0.00; y_offset =-0.001 # this works for #'kappa' as y
# #     x_offset = 0.0001; y_offset =+0.0005 # this works for #'Delta_Propensity_Old-Young'
# #     ax[0].text(df_data[var_x][line]+x_offset, df_data[var_y][line]+y_offset, \
# #                df_data['Human'][line], horizontalalignment='left', size='small', color='black', weight='regular')
# # get the separating hyperplane in original coordinate
# w = clf.coef_[0]
# a = -w[0]/w[1]
# xx = X[:,0]
# yy = a*xx - (clf.intercept_[0])/w[1]
# raw_plane = scal.inverse_transform([(i,j) for i, j in zip(xx,yy)])
# # plot the separating hyperplane/lines
# ax[0].plot(raw_plane[:,0], raw_plane[:,1], 'k-')
# # ax[0].text(0.05,0.9, 'accuracy = %.2f%%\ncoverage = %.2f%%'%(accuracy_all*100,frac_mapped*100), transform=ax.transAxes)
# # add legend etc.
# # ax.legend(markerscale=6, bbox_to_anchor=(1.4, 1))
# plt.tight_layout()
# print('Accuracy of SVM classifier with %s and %s : %.2f'%(var_x, var_y,accuracy_all))
# #### plot the population
# pop_dir = 'Properties'
# df_pop = pd.read_csv(os.path.join(pop_dir, 'kf_combined_prop_metrics.csv'))
# df_AGG = df_pop[df_pop['Sample.Type']=='AGG']
# df_TL = df_pop[df_pop['Sample.Type']=='TL']
# key = ['N. furzeri Protein Id','Tissue']
# df_pos_AGG = df_AGG[(df_AGG['OvY_logFC']>0) & (df_AGG['OvY_pval']<=0.05)]
# df_pos_TL = df_TL[(df_TL['OvY_logFC']>0) & (df_TL['OvY_pval']<=0.05)]
# pos_AGG_key = df_pos_AGG['Tissue']+df_pos_AGG['N. furzeri Protein Id']
# pos_TL_key = df_pos_TL['Tissue']+df_pos_TL['N. furzeri Protein Id']
# df_pos_AGG['pos_TL'] = pos_AGG_key.isin(pos_TL_key) 
# # plot original data
# sns.scatterplot(x = var_x, y = var_y, style = 'Tissue', hue='pos_TL',\
#                 data = df_pos_AGG ,ax = ax[1],size='OvY_logFC',sizes=(0.1, 500))
# ax[1].plot(raw_plane[:,0], raw_plane[:,1], 'k-')
# fig.savefig(os.path.join(out_dir,'Pop_svm.pdf'))
# #################### End of Section 6 ###################

# # #################### Section 7: Visualize the parameter space I explored ####################
# sns.set(font_scale = 2)
# sns.set_style("white")
# f_dir = 'SysSeeding/SVM'
# df_svm = pd.read_csv(os.path.join(f_dir, 'SysSeedingLinearSVC.csv'))
# df_svm.rename(columns={'Score_100mean':'Score_mean'}, inplace=True)
# plot_coln = ['var_x','var_y','Score_mean']
# df_svm2 = df_svm[plot_coln].copy()
# df_svm2.columns=['var_y','var_x','Score_mean']
# df_plot = df_svm[plot_coln].append(df_svm2, ignore_index = True)
# print (df_plot.shape)
# df_plot.drop_duplicates(inplace=True)
# print (df_plot.shape)
# # fig_svm, ax_svm = plt.subplots(1,1, figsize=(12,12))
# df_plot = df_plot.pivot(*plot_coln)
# df_plot.fillna(0, inplace=True)
# fig_svm = sns.clustermap(df_plot, vmax=100, square=True, figsize=(16, 16))
# fig_svm.savefig(os.path.join(f_dir, 'SVM_score.pdf'))
# #################### End of Section 7 ###################

#################### Section 8: Pairwise correlation on all metrics ####################
# seeding_dir = 'SysSeeding'
# df_data = pd.read_csv(os.path.join(seeding_dir, 'SysSeedingAll_metrics.csv'))
# features = ['OvY_logFC', 'OvY_prop_logFC','log2_Young_TL',\
#             'AA Length', 'FImaxrun','Longest Disorder Stretch',\
#             'FracExpanding','FracDisoPromoting', 'Disordered Frac','Frac_FInumaa',\
#             'PAPAprop','NLLR','MWScore','Frac_QN',\
#             'pI', 'FracNeg','Frac_pos','FracCharged',\
#             'delta', 'kappa','Omega','MeanNetCharge','uversky_hydropathy','Frac_neu',\
#             'Frac_aromatic','Frac_aliphatic']
# pmatrix_out = df_data[features].corr(method = 'pearson')
# pmatrix_out.to_csv(os.path.join(seeding_dir, 'SysSeedingLinear_r.csv'))
# fig_r = sns.clustermap(pmatrix_out, vmax=1, square=True, center=0, figsize=(16, 16))
# fig_r.savefig(os.path.join(seeding_dir, 'SysSeedingLinear_r.pdf'))
#################### End of Section 8 ###################
