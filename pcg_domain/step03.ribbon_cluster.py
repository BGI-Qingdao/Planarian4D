#!/usr/bin/env python3
###############################################################################
### import package
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
import umap.plot
from scipy.signal import savgol_filter
from sklearn.mixture import BayesianGaussianMixture
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MaxAbsScaler
from sklearn.svm import SVC
from sklearn import tree
import seaborn as sns
from sklearn.linear_model import LogisticRegression
###############################################################################
### load raw data
pcgs = np.loadtxt('names',dtype=str)

samples = [ 'WT', '0hpa1', '0hpa2',  '12hpa1', '12hpa2', '36hpa1', '36hpa2', 
            '3dpa1', '3dpa2' ,'5dpa1', '5dpa2', '7dpa1', '7dpa2',
            '10dpa1', '10dpa2' , '14dpa1', '14dpa2' ]

raw_data = pd.read_csv('scaled_smooth_sigmod2.csv',sep=',',header=0)
raw_data.columns = ['idx','sample','AP','t_min','t_max','exp_cell_count','exp_sum','cell_count','density1k','gene','rtype']
raw_data.drop(columns = ['idx'],inplace=True)
sample_data = {}
for sp in samples:
    sample_data[sp] = raw_data[ raw_data ['sample'] == sp ].copy().reset_index()

###############################################################################
### construct X matrix
sample_header = {} 
sample_X = {}
for sp in samples:
    sample_raw_data = sample_data[sp]
    header = None 
    features = []
    for smes in pcgs:
        if header is None:
             header = sample_raw_data[sample_raw_data['gene'] == smes]
             header = header[['sample','AP','rtype']].copy()
             header = header.reset_index()
             sample_header[sp] = header
        raw_v = sample_raw_data[sample_raw_data['gene'] == smes]['density1k'].copy()
        features.append(raw_v)
    all_features = np.vstack(features)
    sample_X[sp] = all_features.T

#temp_X_list = []
#for sp in samples:
#     temp_X_list.append(sample_X[sp])
#temp_X = np.concatenate(temp_X_list, axis=0)
sample_scaled_X = sample_X

###############################################################################
### do PCA
sample_scaled_pca_X = {}
pca =  PCA()
pca.fit(sample_scaled_X['WT'])
for sp in samples:
    sample_scaled_pca_X[sp] = pca.transform(sample_scaled_X[sp])
    sample_header[sp]['pc0'] = sample_scaled_pca_X[sp][:,0]
    sample_header[sp]['pc1'] = sample_scaled_pca_X[sp][:,1]
### draw pca analysis
pca_component_weight = pca.explained_variance_ratio_
pca_index = np.arange(1,len(pca_component_weight)+1)
plt.figure(figsize=(8,6))
plt.plot(pca_index,pca_component_weight)
plt.xlabel('PCs')
plt.ylabel('variance ratio')
plt.savefig('pca_component_WT.png',dpi=100)
plt.close()

pca_component_weight = np.cumsum(pca_component_weight)
plt.figure(figsize=(8,6))
plt.plot(pca_index,pca_component_weight)
plt.xlabel('PCs')
plt.ylabel('Cumulative variance ratio')
plt.savefig('pca_component_cumsum_WT.png',dpi=100)
plt.close()

plt.figure(figsize=(8,6))
plt.scatter(sample_scaled_pca_X['WT'][:,0] ,sample_scaled_pca_X['WT'][:,1])
plt.xlabel('PC0')
plt.ylabel('PC1')
plt.savefig('pca_space_WT.png',dpi=100)
plt.close()

###############################################################################
### do UMAP
sample_scaled_umap_X = {}
umapfit = umap.UMAP(metric='manhattan',random_state=0)
umap_mapper = umapfit.fit(sample_scaled_X['WT'])
for sp in samples:
    sample_scaled_umap_X[sp] = umap_mapper.transform(sample_scaled_X[sp])
    sample_header[sp]['umap0'] = sample_scaled_umap_X[sp][:,0]
    sample_header[sp]['umap1'] = sample_scaled_umap_X[sp][:,1]
### draw umap analysis
umap.plot.points(umap_mapper)
plt.savefig(f"umap_space_WT.png")

###############################################################################
## do GM in WT and then label transfer to others by SVM
for n_cluster in [ 3, 4,5,6,7,8,9,10,11,12,13,14,15 ] :
    print(f'n_cluster={n_cluster}')
    clusterer = BayesianGaussianMixture(n_components=n_cluster, random_state=0,covariance_type="diag")
    labels = clusterer.fit_predict(sample_scaled_X['WT'])
    sample_header['WT']['label'] = labels
    #DTclf = tree.DecisionTreeClassifier()
    #DTclf = DTclf.fit(sample_scaled_X['WT'],labels)
    #svc = SVC()
    lr = LogisticRegression()
    lr.fit(sample_scaled_X['WT'],labels)
    for sp in samples:
        if sp == 'WT':
            continue
        labels = lr.predict(sample_scaled_X[sp])
        #labels = DTclf.predict(sample_scaled_X[sp])
        sample_header[sp]['label'] = labels
    results = []
    for sp in samples:
        results.append(sample_header[sp])
    current_result = pd.concat(results, ignore_index=True)
    current_result.to_csv(f'lt_lr_{n_cluster}.csv',index=None,sep=',')

