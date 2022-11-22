#!/usr/bin/env  python3

import sys
import pandas as pd
import numpy as np
from skimage import io as skio

def draw_sample( sample, sample_info ,binsize, color_list, draw_binsize = 10):
    step = int( binsize / draw_binsize )
    pos_file = f'cell_pos_db/{sample}.txt'
    pos_data = pd.read_csv( pos_file, sep=' ', header=None )
    pos_data.columns = ['x','y','z']
    pos_data = pos_data.groupby(['x', 'y']).agg(z=('z', 'max')).reset_index()
    
    xmin = np.min( pos_data['x'] )
    ymin = np.min( pos_data['y'] )
    pos_data['x'] = pos_data['x'] - xmin 
    pos_data['y'] = pos_data['y'] - ymin 
    pos_data['x'] = pos_data['x'] / draw_binsize
    pos_data['y'] = pos_data['y'] / draw_binsize
    pos_data = pos_data.astype(int)
    #y,x = np.nonzero(pos_data)
    W = np.max( pos_data['y'] ) + 10
    H = np.max( pos_data['x'] ) + step 

    H_min = int((-100 - ymin )/binsize)
    H_max = int((100 - ymin )/binsize)+1

    drawing_data = np.zeros((H,W,3),dtype='uint8')
    body_data = np.zeros((H,W,3),dtype='uint8')
    body_data[ pos_data['x'],pos_data['y'] , : ] = 125
    ret_data = np.zeros((H,W,3),dtype='uint8')
    ################
    # draw mid
    sample_label = sample_info[ sample_info['rtype'] == 'mid' ]['label'].to_numpy()
    sample_label = sample_label.astype(int)
    ap_list = sample_info[ sample_info['rtype'] == 'mid' ]['AP'].to_numpy()
    ap_list_new = ap_list - int(xmin / draw_binsize)
    for j,i in enumerate(ap_list_new):  #range(len(sample_label)):
        drawing_data[i*step:(i+1)*step,:,:] = 125
        if sample_label[j] != -1 :
            drawing_data[i*step:(i+1)*step,H_min:H_max,0] = color_list[sample_label[j],0]
            drawing_data[i*step:(i+1)*step,H_min:H_max,1] = color_list[sample_label[j],1]
            drawing_data[i*step:(i+1)*step,H_min:H_max,2] = color_list[sample_label[j],2]
        else:
            drawing_data[i*step:(i+1)*step,:,:] = 50
            
        body_data[ i*step:(i+1)*step,:,: ] = 0
    ################
    # draw left
    sample_label = sample_info[ sample_info['rtype'] == 'other' ]['label'].to_numpy()
    sample_label = sample_label.astype(int)
    ap_list = sample_info[ sample_info['rtype'] == 'other' ]['AP'].to_numpy()
    ap_list_new = ap_list - int(xmin / draw_binsize)
    for j,i in enumerate(ap_list_new):  #range(len(sample_label)):
        if sample_label[j] != -1 :
            drawing_data[i*step:(i+1)*step,:H_min,0] = color_list[sample_label[j],0]
            drawing_data[i*step:(i+1)*step,:H_min,1] = color_list[sample_label[j],1]
            drawing_data[i*step:(i+1)*step,:H_min,2] = color_list[sample_label[j],2]

    ################
    # draw right
    sample_label = sample_info[ sample_info['rtype'] == 'other' ]['label'].to_numpy()
    sample_label = sample_label.astype(int)
    ap_list = sample_info[ sample_info['rtype'] == 'other' ]['AP'].to_numpy()
    ap_list_new = ap_list - int(xmin / draw_binsize)
    for j,i in enumerate(ap_list_new):  #range(len(sample_label)):
        if sample_label[j] != -1 :
            drawing_data[i*step:(i+1)*step,H_max:,0] = color_list[sample_label[j],0]
            drawing_data[i*step:(i+1)*step,H_max:,1] = color_list[sample_label[j],1]
            drawing_data[i*step:(i+1)*step,H_max:,2] = color_list[sample_label[j],2]
    
    ret_data[ pos_data['x'],pos_data['y'] , : ] = drawing_data[  pos_data['x'],pos_data['y'] , : ]    
    ret_data[ np.nonzero(body_data) ] = body_data[ np.nonzero(body_data) ]
    skio.imsave(f'{prefix}.{sample}.tif',ret_data)
    return ret_data

filename = sys.argv[1]
binsize = int(sys.argv[2])
prefix = sys.argv[3]

color_list = np.array([ 
                [255,0  ,0  ],  # type 0 <--> red
                [  0,255,0  ],  # type 1 <--> green
                #[255,  0,255],  # type 1 <--> cyan
                [255,255,0  ],  # type 2 <--> yellow
                [  0,255,255],  # type 3 <--> cyan 
                [255,  0,255],  # type 4 <--> magenta 
                [255,255,255],  # type 5 <--> white
                [  0,128,255],  #type 6 <--> blue-like
                [255,246,143],  # type 7 <--> Khaki1 
                [255,231,186],  # type 8 <--> Wheat
                [255,20 ,147],  # type 9 <--> DeepPink
                [205,133,63 ],  # type 10 <-> Peru	
                [250,128,114],  # type 11 <-> Salmon
                [224,238,238],  # type 12 <-> Azure2
             ])

cluster_result = pd.read_csv(filename, sep=',',header=0)
#print(cluster_result.head(),flush=True)
H_max = 0
W_total = 0
imgs = []
for sample in [ 'WT', '0hpa1', '0hpa2', '12hpa1', '12hpa2', '36hpa1', '36hpa2', '3dpa1', '3dpa2', '5dpa1', '5dpa2', '7dpa1', '7dpa2', '10dpa1', '10dpa2','14dpa1', '14dpa2' ] :
    sample_info= (cluster_result[cluster_result['sample'] == sample][['AP','rtype','label']])
    draw_img = draw_sample( sample, sample_info, binsize, color_list )
    H , W , _= draw_img.shape
    W_total += W
    imgs.append(draw_img)
    if H > H_max :
        H_max = H

final_img = np.zeros((H_max,W_total,3),dtype='uint8')
W_start = 0
for img in imgs:
    H , W ,_ = img.shape
    final_img[:H,W_start:W_start+W,:] = img[:,:,:]
    W_start = W_start + W

skio.imsave(f'{prefix}.png',final_img)
