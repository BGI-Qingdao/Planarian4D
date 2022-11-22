#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from matplotlib import rcParams
import seaborn as sns
from scipy.ndimage import gaussian_filter1d

###########################################################################################
## sigmoid
def sigmoid(x):
  x1 = x.copy()
  x = x - 2 
  return 1 / (1 + np.exp(-x))

###########################################################################################
## load data
gene_names = pd.read_csv('names',sep='\t',header=None)
gene_names.columns=['SMESID','name']
sample_names =  ['WT','0hpa1','0hpa2','12hpa1','12hpa2','36hpa1','36hpa2','3dpa1','3dpa2',
                '5dpa1','5dpa2','7dpa1','7dpa2','10dpa1','10dpa2','14dpa1','14dpa2' ]

raw_data_mid = pd.read_csv('raw_all.ml100.csv',sep=',',header=0, index_col=0)
raw_data_other = pd.read_csv('raw_all.ml100.other.csv',sep=',',header=0, index_col=0)

###########################################################################################
## prepare data
raw_data_mid['rtype'] = ['mid'] * len(raw_data_mid)
raw_data_other['rtype'] = ['other'] * len(raw_data_other)
###########################################################################################
## smooth
for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data_mid[ raw_data_mid['gene'] ==  smedid]
     for x in sample_names:
         raw = drawdata[drawdata['sample'] == x ]['density1k'].to_numpy()
         raw = gaussian_filter1d(raw,sigma=3)
         raw_data_mid.loc[((raw_data_mid['gene'] ==  smedid)&(raw_data_mid['sample'] == x)), 'density1k'] = raw


for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data_other[ raw_data_other['gene'] ==  smedid]
     for x in sample_names:
         raw = drawdata[drawdata['sample'] == x ]['density1k'].to_numpy()
         raw = gaussian_filter1d(raw,sigma=3)
         raw_data_other.loc[((raw_data_other['gene'] ==  smedid)
                           &(raw_data_other['sample'] == x)), 'density1k'] = raw

raw_data = pd.concat((raw_data_mid,raw_data_other),ignore_index = True)
raw_data.to_csv('smooth.csv')

for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data[ raw_data['gene'] ==  smedid ]
     sns.relplot(data=drawdata,x='AP',y='density1k',style='rtype', col='sample', col_wrap=6, kind="line",)
     plt.title(name)
     plt.savefig(f'{name}.smooth.png')
     plt.close()
###########################################################################################
## scale and center
for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata0 = raw_data_mid[ raw_data_mid['gene'] ==  smedid]
     drawdata1 = raw_data_other[ raw_data_other['gene'] ==  smedid]
     
     raw0 = drawdata0['density1k'].to_numpy()
     raw1 = drawdata1['density1k'].to_numpy()
     raw = np.hstack((raw0,raw1))
     mean = np.mean(raw)
     std = np.std(raw)
     raw = raw - mean
     raw = raw / std
     
     l1 = len(raw0)
     l2 = len(raw1)
     raw0 = raw[:l1]
     raw1 = raw[l1:]
     raw_data_mid.loc[raw_data_mid['gene'] ==  smedid, 'density1k'] = raw0
     raw_data_other.loc[raw_data_other['gene'] ==  smedid, 'density1k'] = raw1

raw_data = pd.concat((raw_data_mid,raw_data_other),ignore_index = True)
raw_data.to_csv('scaled.smoothed.csv')

for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data[ raw_data['gene'] ==  smedid]
     sns.relplot(data=drawdata,x='AP',y='density1k',style='rtype', col='sample', col_wrap=6, kind="line",)
     plt.title(name)
     plt.savefig(f'{name}.scaled_smoothed.png')
     plt.close()



###########################################################################################
## sigmoid
for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data[ raw_data['gene'] ==  smedid]
     raw = drawdata['density1k'].to_numpy()
     raw = sigmoid(raw)
     raw_data.loc[raw_data['gene'] ==  smedid, 'density1k'] = raw

for index , row in gene_names.iterrows():
     name = row['name']
     smedid = row['SMESID']
     drawdata = raw_data[ raw_data['gene'] ==  smedid]
     sns.relplot(data=drawdata,x='AP',y='density1k', style='rtype', col='sample', col_wrap=6, kind="line",)
     plt.title(name)
     plt.savefig(f'{name}.scaled.smooth.sigmod2.png')
     plt.close()
raw_data.to_csv('scaled_smooth_sigmod2.csv')
