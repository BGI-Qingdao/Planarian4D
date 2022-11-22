#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

def count_by_center(raw_data , xid,yid, binsize):
    return np.sum ( (raw_data.x >= (xid - binsize)) & (raw_data.x <=  (xid +binsize)) & (raw_data.y >=(yid - binsize)) & (raw_data.y <=(yid + binsize)) )

####################################################
# args
infile = sys.argv[1]
prefix = sys.argv[2]
binsize = 3 #int(sys.argv[3])  # binsize in 5 micron =>15 micron

print(f'expand rect_edge = {binsize}*2 ',flush=True)
####################################################
# load raw data
#raw_data = pd.read_csv(infile,sep='\t',header=None)
#raw_data.columns = ['slice','cell','x','y','z','c']
raw_data = pd.read_csv(infile,skipinitialspace=True,header = None, delim_whitespace=True)
raw_data.columns = ['x','y','z','c','c2']

raw_data['x'] = raw_data['x'] / 5
raw_data['y'] = raw_data['y'] / 5 

raw_data['x'] = raw_data['x'].astype(int)
raw_data['y'] = raw_data['y'].astype(int)
raw_data['c'] = raw_data['c'].astype(int)
raw_data['c2'] = raw_data['c2'].astype(int)

###################################################
# get 2D idx for calc density
xmin = np.min(raw_data['x'])
xmax = np.max(raw_data['x'])
ymin = np.min(raw_data['y'])
ymax = np.max(raw_data['y'])

width = xmax - xmin +1
height = ymax - ymin +1

raw_data['x'] = raw_data['x'] - xmin
raw_data['y'] = raw_data['y'] - ymin

mask = np.zeros((height,width),dtype=int)
mask[ raw_data['y'], raw_data['x'] ] = 1
yidx ,xidx = mask.nonzero()

###################################################
# get total number of cells
total = np.zeros((height,width),dtype=float)
sumcell = 0
for i in range(len(yidx)):
    yid, xid = yidx[i] ,xidx[i]
    total[ yid, xid ] = count_by_center(raw_data , xid,yid, binsize)
    sumcell += total[ yid, xid ]

print(f'total cell : {sumcell}',flush=True)

for j in range(1,10,1):
    sumcell = 0
    ct_data = raw_data[raw_data['c'] == j].copy()
    count = np.zeros((height,width),dtype=float)
    for i in range(len(yidx)):
        yid, xid = yidx[i] ,xidx[i]
        count[ yid, xid ] =  count_by_center( ct_data,xid ,yid, binsize)
        sumcell += count[ yid, xid ]
    print(f'total cell of {j}: {sumcell}',flush=True)
    for i in range(len(yidx)): 
        yid, xid = yidx[i] ,xidx[i]
        count[ yid, xid ] = count[ yid, xid ] / total[ yid, xid ]
    tsd_count = count * 10000
    tsd_count = tsd_count.astype(int)
    y_coords , x_coords = np.nonzero(tsd_count)
    density_df = pd.DataFrame()
    density_df['x'] = x_coords
    density_df['y'] = y_coords
    density_df['z'] = y_coords #fake z
    density_df['x'] = density_df['x'] + xmin
    density_df['x'] = density_df['x'] * 5
    density_df['y'] = density_df['y'] + ymin
    density_df['y'] = density_df['y'] * 5
    density_df['c'] = tsd_count[y_coords,x_coords]

    density_df.to_csv(f'{prefix}.{j}.density_list.txt',header=None,index=None,sep=' ')
    #np.savetxt(f'{prefix}.{j}.density.txt',count)

