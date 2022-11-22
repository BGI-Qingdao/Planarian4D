#!/usr/bin/env  python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('dark_background')

pcgs = np.loadtxt('name_map',delimiter='\t',dtype=str)

raw_data = pd.read_csv('scaled_smooth_sigmod2.csv',sep=',',header=0, index_col=0)
ret_data = pd.read_csv('lt_lr_10.csv',sep=',',header=0)

len_array = len(pcgs)

target_df = ret_data

for i in range(len_array):
    smes = pcgs[i,1] 
    name = pcgs[i,0] 
    one_gene_data = raw_data[ raw_data['gene'] == smes ].copy()
    target_df[name] = one_gene_data['density1k'].to_numpy()

target_df.to_csv('lt_lr_10.flat_genes.csv')
target_df['label'] = target_df['label'].astype(str)
target_df = target_df[target_df['rtype']=="mid"].copy()
target_df = target_df[target_df['sample']=="WT"].copy()
target_df = target_df.drop(columns=['index','sample','AP','pc0','pc1','umap0','umap1','rtype'])
target_df = target_df.groupby(['label']).mean()
target_df.to_csv('heatmap.csv')

AP_list = []
AP_list.append('H-polar1')
AP_list.append('H-polar2')
AP_list.append('Neck')
AP_list.append('Trunk')
AP_list.append('P-start')
AP_list.append('P-mid')
AP_list.append('P-end')
AP_list.append('Tail1')
AP_list.append('Tail2')
AP_list.append('T-polar')

names = {}
names[5] = 'H-polar1'
names[7] = 'H-polar2' 
names[3] = 'Neck'
names[1] = 'Trunk'
names[4] = 'P-start'
names[8] = 'P-mid'
names[0] = 'P-end'
names[9] = 'Tail1'
names[6] = 'Tail2'
names[2] = 'T-polar'

newhead = []
for i in range(10):
    newhead.append(names[i])

#newhead = [ 'Pharynx end' ,'Trunk','Head','Pharynx middle','Tail polar','Pharynx start','Tail' ]

print(target_df.index,flush=True)
target_df['new_names'] =newhead 
target_df = target_df.set_index(['new_names'])
#newhead = ['Head','Trunk','Pharynx start','Pharynx middle', 'Pharynx end' ,'Tail','Tail polar']
target_df = target_df.reindex(AP_list)

target_df = target_df[[ 'sfrp-1','fz5_8-3','ndl-4','ndl-1','fz5_8-4','ndk','zic-2','prep','wnt2','ndl-5','sfrp-2','ndl-2','ndl-3','wntP-2_wnt11-5','ptk7','teashirt','wntA','dd4427','dd13065','sp5','fz4-1','Post-2c','hox4b','wnt11-2','wntless']].copy()
sns.heatmap(target_df, linewidths=2,cmap='Spectral_r',vmin=0,vmax=1)
plt.show()
