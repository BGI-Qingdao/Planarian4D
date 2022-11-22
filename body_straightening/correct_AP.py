#!/usr/bin/env python3
import math
import os
from skimage import io
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.spatial import KDTree
import argparse

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
    c = pd.read_csv('constrain.txt', sep='\t')
    c = c.set_index('id')
    p = pd.read_csv('line_parameter.txt', sep='\t')
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
    #raw.columns = ['id','cell','x','y','z']
    ## get up and down
    #raw = pd.DataFrame()
    #raw['x'] = x
    #raw['y'] = y
    raw0 = raw[(raw['x']<=x2)&(raw['x']>=x1)].copy()
    raw1 = raw[raw['x']<x1].copy()
    raw2 = raw[raw['x']>x2].copy()
    raw0['y_inline'] = np.polyval(k, raw0['x'])
    raw1['y_inline'] = np.polyval(k1, raw1['x'])
    raw2['y_inline'] = np.polyval(k2, raw2['x'])
    raw = pd.concat([raw0,raw1,raw2], axis=0, ignore_index=True)
    raw['up'] = raw['y']>=raw['y_inline']
    raw['down'] = raw['y']<raw['y_inline']

    # get clost point for all raw points
    kdtree = KDTree(np.array(line_new_points))
    raw_pos = np.zeros((len(raw),2))
    raw_pos[:,0] = raw['x']
    raw_pos[:,1] = raw['y']
    darray, iarray = kdtree.query(raw_pos)

    # get new nx ny
    raw['ny'] = darray
    new_y = np.zeros(len(raw['ny']))
    new_y[raw['up'] ]  = raw[raw['up']]['ny']
    new_y[raw['down']] = -raw[raw['down']]['ny']
    raw['ny'] = new_y

    raw['nx'] = length[np.array(iarray).astype(int)]

    raw['nx'] = raw['nx'].astype(int)
    raw['ny'] = raw['ny'].astype(int)
    #raw.to_csv(f'{args.output}/index_pos_{name}.csv',sep=',')
    new = pd.DataFrame()
    new[['id','cell','x','y','z']] = raw[['id','cell','nx','ny','z']]
    '''if not os.path.exists(args.output):
        os.mkdir(args.output)
    new.to_csv(f'{args.output}/correct_pos_{name}.csv',sep=',',header=False, index=False,index_label=False)'''
    # draw final
    '''raw['z'] = raw['z']/10
    raw['z'] = raw['z']-raw['z'].min()
    raw['z'] = raw['z'].astype(int)
    min_y = raw['ny'].min()
    raw['nx'] = raw['nx']-raw['nx'].min()
    raw['ny'] = raw['ny']-raw['ny'].min()
    
    data = np.zeros((raw['z'].max()+1,raw['ny'].max()+10,raw['nx'].max()+10))
    data[raw['z'],raw['ny'],raw['nx']] =255
    data[raw['z'],-min_y:-min_y+5,:]=255
    data[raw['z'],-min_y-5:-min_y,:]=255
    data = data.astype('uint8')

    io.imsave(f'{args.output}/correct_{name}.tif',data)'''

    #rotation
    #coor = pd.read_csv(f'{args.output}/correct_pos_{name}.csv', sep=',', names=['id', 'cell', 'x', 'y', 'z'])
    coor = new.copy()
    r = pd.read_csv("/dellfsqd2/ST_OCEAN/USER/huangzhi/wt/correct_AP/new_pos/rotation.txt", sep='\t', names=['id', 'rotation'])
    r = r.set_index('id')
    rotation = r.loc[name, 'rotation']
    if rotation == 'TRUE':
        coor['x'] = -coor['x']
        shift = -(coor[coor['y'] == 0]['x'].min())
        coor['x'] = coor['x'] + shift
        coor['y'] = -coor['y']

    else:
        shift = -(coor[coor['y'] == 0]['x'].min())
        coor['x'] = coor['x'] + shift
        # coor['y'] = -coor['y']

    if name == '10dpa1' or name == '36hpa2':
        z_max = coor['z'].max()
        z_min = coor['z'].min()
        coor['z'] = -coor['z']
        coor['y'] = -coor['y']
        coor['z'] = coor['z'] + z_max + z_min
        #Zshift_list = Zshift_list.append({'id': name, 'Zshift': z_max + z_min}, ignore_index=True)

    #P = coor[coor['y'] == 0]['x'].max()
    #P_coor = P_coor.append({'id': name, 'P_x': P}, ignore_index=True)
    #shift_list = shift_list.append({'id': name, 'shift': shift}, ignore_index=True)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    coor.to_csv(f'{args.output}/correct_pos_{name}.csv', sep=',', header=False, index=False, index_label=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", help="the name of sample,like WT", type=str, default='', nargs='?')
    parser.add_argument("-p", "--posfile", help="input pos file path", type=str, default='', nargs='?')
    #parser.add_argument("-c", "--constrain", help="constrain file path", type=str, default='', nargs='?')
    #parser.add_argument("-d", "--define", help="Domain definition, x_max x_min y_max y_min", type=str, default=[], nargs='+')
    parser.add_argument("-o", "--output", help="directory to save files", type=str,default='APoutput')
    args = parser.parse_args()
    main(args)
