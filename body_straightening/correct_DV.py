#!/usr/bin/env python3
import math
from skimage import io
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.spatial import KDTree
import argparse
import os

def getCurveY(xvals0,xvals1,xvals2,x1,x2,y1,y2,k):
    yvals0 = np.polyval(k, xvals0)

    k1 = k[2]+2*k[1]*x1+3*k[0]*x1*x1
    b1 = y1 - k1*x1
    k1 = [k1,b1]
    yvals1 = np.polyval(k1, xvals1)

    k2 = k[2]+2*k[1]*x2+3*k[0]*x2*x2
    b2 = y2 - k2*x2
    k2 = [k2,b2]
    yvals2 = np.polyval(k2, xvals2)

    yvals = np.hstack((yvals1, yvals0, yvals2))

    return yvals,k1,k2

def main(args):
    #load parameter and constrain
    c = pd.read_csv("constrain.txt", sep='\t')
    c = c.set_index('id')
    p = pd.read_csv("line_parameter.txt", sep='\t')
    p = p.set_index('id')
    name = args.name
    #load line
    x1 = c.loc[name, 'x1']*5
    x2 = c.loc[name, 'x2']*5
    y1 = c.loc[name, 'y1']*5
    y2 = c.loc[name, 'y2']*5
    k = [p.loc[name, 'k0']*0.2*0.2, p.loc[name, 'k1']*0.2, p.loc[name, 'k2'], p.loc[name, 'k3']*5]

    # line_new_points
    xvals0 = np.arange(x1,x2,0.1)
    xvals1 = np.arange(0,x1,0.1)
    xvals2 = np.arange(x2,x2+500,0.1)
    x = np.hstack((xvals1,xvals0,xvals2))

    line_new_points = np.zeros((len(x),2))
    line_new_points[:,0] = x
    y,k1,k2 = getCurveY(xvals0,xvals1,xvals2,x1,x2,y1,y2,k)
    line_new_points[:,1] = y

    # calc distance of curve as new x
    length = np.zeros(len(line_new_points))
    for i in range(1,len(line_new_points)):
        length[i] = distance.euclidean(line_new_points[i,:],line_new_points[i-1,:])
    length = np.cumsum(length)

    # load raw points
    raw = pd.read_csv(args.posfile,sep=",",names=['id','cell','x','y','z'])
    raw['x'] = raw['x'] + 150
    #raw.columns = ['id','cell','x','y','z']
    ## get up and down
    #raw = pd.DataFrame()
    #raw['x'] = x
    #raw['y'] = y
    raw0 = raw[(raw['x']<=x2)&(raw['x']>=x1)].copy()
    raw1 = raw[raw['x']<x1].copy()
    raw2 = raw[raw['x']>x2].copy()
    raw0['z_inline'] = np.polyval(k, raw0['x'])
    raw1['z_inline'] = np.polyval(k1, raw1['x'])
    raw2['z_inline'] = np.polyval(k2, raw2['x'])
    raw = pd.concat([raw0,raw1,raw2], axis=0, ignore_index=True)
    raw['up'] = raw['z']>=raw['z_inline']
    raw['down'] = raw['z']<raw['z_inline']

    # get clost point for all raw points
    kdtree = KDTree(np.array(line_new_points))
    raw_pos = np.zeros((len(raw),2))
    raw_pos[:,0] = raw['x']
    raw_pos[:,1] = raw['z']
    darray, iarray = kdtree.query(raw_pos)

    # get new nx ny
    raw['nz'] = darray
    new_z = np.zeros(len(raw['nz']))
    new_z[raw['up'] ]  = raw[raw['up']]['nz']
    new_z[raw['down']] = -raw[raw['down']]['nz']
    raw['nz'] = new_z

    raw['nx'] = length[np.array(iarray).astype(int)]

    raw['nx'] = raw['nx'].astype(int)
    raw['nz'] = raw['nz'].astype(int)
    #raw.to_csv(f'{args.output}/index_pos_{name}.csv',sep=',')
    new = pd.DataFrame()
    new[['id','cell','x','y','z']] = raw[['id','cell','nx','y','nz']]
    x_shift = new[new['y']==0]['x'].min()
    new['x'] = new['x'] - x_shift
    print(name, x_shift)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    new.to_csv(f'{args.output}/correct_pos_{name}.csv',sep=',',header=False, index=False,index_label=False)
    # draw final
    '''raw['nz'] = raw['nz']/10
    raw['nz'] = raw['nz'].astype(int)
    min_z = raw['nz'].min()
    raw['nz'] = raw['nz']-raw['nz'].min()
    raw['nx'] = raw['nx']-raw['nx'].min()
    raw['y'] = raw['y']-raw['y'].min()
    
    data = np.zeros((raw['nz'].max()+1,raw['y'].max()+10,raw['nx'].max()+10))
    data[raw['nz'],raw['y'],raw['nx']] =255
    data[-min_z,:,:]=255
    data[-min_z,:,:]=255
    data = data.astype('uint8')
    io.imsave(f'{args.output}/correct_{name}.tif',data)'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", help="the name of sample,like WT", type=str, default='', nargs='?')
    parser.add_argument("-p", "--posfile", help="input pos file path", type=str, default='', nargs='?')
    #parser.add_argument("-c", "--constrain", help="constrain file path", type=str, default='', nargs='?')
    #parser.add_argument("-d", "--define", help="Domain definition, x_max x_min y_max y_min", type=str, default=[], nargs='+')
    parser.add_argument("-o", "--output", help="directory to save files", type=str,default='DVoutput')
    args = parser.parse_args()
    main(args)
