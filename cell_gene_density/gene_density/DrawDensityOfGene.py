#!/usr/bin/env python3

import sys
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#############################################################
# get parameters
#############################################################
gene_smes = sys.argv[1]
if len(sys.argv) > 2 :
    binsize = int(sys.argv[2])
else :
    binsize = 20

#############################################################
# load data
#############################################################
WT_pos = '/node02/guolidong/wochong/DensityPlot/WT.txt'
WT_pos_data = pd.read_csv(WT_pos, sep=' ',header=None)
WT_pos_data.columns = ['x','y','z']

gene_expr_file = f'/node02/guolidong/wochong/DensityPlot/gene_atlas_smes_all/WT/{gene_smes}.txt'
gene_expr =  pd.read_csv(gene_expr_file,sep=' ',header=None)
gene_expr.columns = ['x','y','z','v']
#############################################################
# construct Body area
#############################################################
xmin = np.min(WT_pos_data['x'])
xmax = np.max(WT_pos_data['x'])
ymin = np.min(WT_pos_data['y'])
ymax = np.max(WT_pos_data['y'])
X, Y = np.mgrid[xmin:xmax:binsize, ymin:ymax:binsize]
positions = np.vstack([X.ravel(), Y.ravel()])

body_mask = np.zeros(X.shape)
WT_pos_data['x'] = WT_pos_data['x'] - xmin
WT_pos_data['x'] = WT_pos_data['x'] / binsize
WT_pos_data['x'] = WT_pos_data['x'].astype(int)
WT_pos_data['y'] = WT_pos_data['y'] - ymin
WT_pos_data['y'] = WT_pos_data['y'] / binsize
WT_pos_data['y'] = WT_pos_data['y'].astype(int)
body_mask[ WT_pos_data['x'],WT_pos_data['y'] ] = 2 

#############################################################
# Gen gene density
#############################################################
values = np.vstack([gene_expr['x'], gene_expr['y']])
kernel = stats.gaussian_kde(values,weights=gene_expr['v'])
Z = np.reshape(kernel(positions).T, X.shape)
Z[body_mask==0]=-1
#############################################################
# Draw density 
#############################################################

fig = plt.figure()
ax = fig.add_subplot(111)
plt.contourf(X, Y, Z, levels=np.linspace(0, Z.max(), 10), cmap=plt.cm.Reds)
plt.contour(X, Y, body_mask, levels=[1], colors="k", linestyles="solid")
plt.xticks([])
plt.yticks([])
ax.set_aspect('equal', adjustable='box')
plt.savefig(f'{gene_smes}.png')
