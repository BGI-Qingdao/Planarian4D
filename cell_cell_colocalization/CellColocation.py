import sys
import numpy as np
import pandas as pd
import igraph as ig
from scipy import stats
import matplotlib.pyplot as plt
import random

# keep result Repeatable
random.seed(42)
###########################################################
# grid functions 
###########################################################

#
# get grid points from one slice ( in sample area ) for probability sample
# 
def grid_in_slice(x,y,binsize):
    pos_data_x = x.copy()
    pos_data_y = y.copy()
    xmin = np.min(pos_data_x) - 2*binsize
    xmax = np.max(pos_data_x) + 2*binsize
    ymin = np.min(pos_data_y) - 2*binsize
    ymax = np.max(pos_data_y) + 2*binsize
    X, Y = np.mgrid[xmin:xmax:binsize, ymin:ymax:binsize]
    body_mask = np.zeros(X.shape)
    pos_data_x =pos_data_x - xmin
    pos_data_x =pos_data_x / binsize
    pos_data_x =pos_data_x.astype(int)
    pos_data_y =pos_data_y - ymin
    pos_data_y =pos_data_y / binsize
    pos_data_y =pos_data_y.astype(int)
    body_mask[ pos_data_x,pos_data_y ] = 1
    X_in_slice  = X[body_mask==1]
    Y_in_slice  = Y[body_mask==1]
    return X_in_slice, Y_in_slice

#
# get grid points from all slices ( in sample area ) for probability sample
# 
def get_3D_grid_sample_points(sample_raw_data):
    positions_for_density  = []
    for sid in np.unique(sample_raw_data['z']):
        slice_data = sample_raw_data[sample_raw_data['z'] == sid]
        x_in_slice = slice_data['x'].to_numpy()
        y_in_slice = slice_data['y'].to_numpy()
        x_in_grid , y_in_grid = grid_in_slice(x_in_slice,y_in_slice,binsize)
        positions_for_density.append(np.vstack([x_in_grid,y_in_grid,(np.ones(x_in_grid.shape)*sid)]).T)
    positions = np.vstack(positions_for_density)
    return positions

###########################################################
# visualise functions 
###########################################################

#
# print result graph in dot format
# generate result by dot -Tpng -o test.png test.dot
#
def print_dot(draw_points,draw_lines,labels,cutoff=0.1):
    ct_colors = ['#d60000', '#e2afaf', '#018700', '#a17569', '#e6a500', '#004b00',
     '#6b004f', '#573b00', '#005659', '#5e7b87', '#0000dd', '#00acc6',
     '#bcb6ff', '#bf03b8', '#645472', '#790000', '#0774d8', '#729a7c',
     '#8287ff', '#ff7ed1', '#8e7b01', '#9e4b00', '#8eba00', '#a57bb8',
     '#5901a3', '#8c3bff', '#a03a52', '#a1c8c8', '#f2007b', '#ff7752',
     '#bac389', '#15e18c', '#60383b', '#546744', '#380000', '#e252ff', ]

    print("strict graph {")
    used_node = []
    for _ ,row in draw_lines.iterrows():
        if row['weight']<cutoff :
            continue
        fromNode = int(row['from'])
        toNode = int(row['to'])
        #cf = row['changefold']
        pw = float(int(row['weight']*100))/10.0
        print(f' {fromNode} -- {toNode} [ penwidth={pw} ]')
        used_node.append(fromNode)
        used_node.append(toNode)
    draw_points = draw_points[draw_points['cid'].isin(used_node)].copy()
    draw_points['fraction'] = draw_points['fraction'].astype(float)
    draw_points['fraction'] = np.log10(draw_points['fraction'])
    draw_points['fraction'] = np.interp( draw_points['fraction'].to_numpy(),(np.min( draw_points['fraction']),np.max( draw_points['fraction'])),(0.25,1.0))
    for _,row in draw_points.iterrows():
        node = int(row['cid'])
        size = row['fraction']
        #print(f' {node} [fillcolor="{color_by_changefold(cf)}" fixedsize=true shape=circle style=filled label={node_label}]')
        print(f' {node} [fillcolor="{ct_colors[node]}" shape=circle style=filled fixedsize=true width={size}]')
    print("}")
    return None

#
# show one 3D density
#
def show_3D_density(x,y,z,density,cid):
    from matplotlib.colors import ListedColormap
    import matplotlib.pylab as pl
    cmap = pl.cm.RdBu
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap[::-1,:])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c=density,cmap=my_cmap,s=1,alpha=0.5)
    w=np.max(x)-np.min(x)
    h=np.max(y)-np.min(y)
    d=np.max(z)-np.min(z)
    ax.set_box_aspect([w,h,d])
    #ax.set_aspect('equal')
    #plt.legend()
    plt.tight_layout()
    plt.savefig(f'3Dtest.{cid}.png')
    plt.close()

###########################################################
# KDE, KL and MST
###########################################################

#
# get sampled density for all grid points
#
def kde3d_all(bootstraped_sample_data,positions, mincell):
    used_ct = []
    Pb_of_allct = []
    for ct in np.unique(bootstraped_sample_data['label']):
        ct_data = bootstraped_sample_data[bootstraped_sample_data['label']==ct]
        cellnum = len(ct_data)
        if cellnum < mincell :
            continue
        used_ct.append(ct)
        all_points = np.vstack([ct_data['x'].to_numpy(),ct_data['y'].to_numpy(),ct_data['z'].to_numpy()])
        kde_ret = stats.gaussian_kde(all_points)
        Pb = kde_ret(positions.T)
        Pb = fill_zero(Pb)
        Pb_of_allct.append(Pb)
    return used_ct,Pb_of_allct
 
#
# fill 0 probability to avoid inf KL value
#
def fill_zero(Pb,esp=1e-20):
    Pb = Pb/np.sum(Pb)
    Pb[Pb<esp] = esp
    return Pb

#
# get KL diversity matrix
# 
def KL(used_ct,Pb_of_allct):  
    kl_array = np.zeros((len(used_ct),len(used_ct)))
    for i,ct1 in enumerate(used_ct):
        for j,ct2 in enumerate(used_ct):
            if ct1 == ct2 :
                continue
            # one row for one ref celltype
            kl_array[i,j] = stats.entropy(pk=Pb_of_allct[i],qk=Pb_of_allct[j],base=2)
    return kl_array

#
# get MST from KL matrix
#
def MSTdf_from_KL(used_ct,kl_array):
    ######################################################
    # get mst
    ######################################################
    yidx,xidx  = np.nonzero(kl_array)
    edges = []
    weights = []
    for yid , xid in zip(yidx,xidx):
        edges.append([yid,xid]) # from KL ref to KL query
        weights.append(kl_array[yid,xid])
    ag = ig.Graph(n=len(used_ct),edges=edges,edge_attrs={'weight':weights},directed=True)
    mst = ag.spanning_tree(weights=ag.es["weight"])
    mst_maskarray = np.array(mst.get_adjacency().data)
    ######################################################
    # gen result
    ######################################################
    yidx,xidx  = np.nonzero(mst_maskarray)
    used_ct = np.array(used_ct,dtype='str')
    round_pd = pd.DataFrame()
    round_pd['from'] = used_ct[yidx.tolist()]
    round_pd['to'] = used_ct[xidx.tolist()]
    return round_pd

###########################################################
# cell type functions
###########################################################
def loading(label_xyz_file,sample):
    raw_meta = pd.read_csv(label_xyz_file,sep='\t',header=0)
    raw_meta['sample'] = raw_meta.apply(lambda x: x['cell'].split('_')[0], axis=1)
    sample_raw_data = raw_meta[raw_meta['sample']==sample].copy()
    sample_raw_data['label'] = sample_raw_data['label'].astype(int)
    total = len(sample_raw_data)
    cids = []
    proportions = []
    for ct in np.unique(sample_raw_data['label']):
        cids.append(ct)
        proportions.append(np.sum(sample_raw_data['label']==ct))
    draw_points = pd.DataFrame()
    draw_points['cid'] = cids
    draw_points['cid']= draw_points['cid'].astype(int)
    draw_points['fraction'] = proportions
    draw_points['fraction'] = draw_points['fraction'].astype(float) 
    return sample_raw_data , draw_points

# main logic
###########################################################
def mainpipe(label_xyz_file,sample,bootnum,binsize):
    # loading meta data
    sample_raw_data, drawpoints = loading(label_xyz_file,sample)
    # create points-in-sample for each slice
    positions = get_3D_grid_sample_points(sample_raw_data)
    # bootstrap
    temp_adjacent_list = []
    for rid in range(bootnum):
        # sample data
        bootstraped_sample_data = sample_raw_data.sample(frac=sample_fac)
        # kde density estimation
        used_ct, Pb_of_allct = kde3d_all(bootstraped_sample_data,positions,mincell)
        if len(used_ct)<2:
           continue
        # KL divergency
        kl_array = KL(used_ct, Pb_of_allct)
        # save KL matrix
        kl_pd = pd.DataFrame(kl_array,index=used_ct,columns=used_ct)
        kl_pd.to_csv(f'{sample}.KL.{rid}.txt',sep='\t')
        # MST
        round_pd = MSTdf_from_KL(used_ct,kl_array)
        round_pd['weight'] = [1.0/float(bootnum)] * len(round_pd)
        temp_adjacent_list.append(round_pd)
    # save result
    final_net = pd.concat(temp_adjacent_list,ignore_index=True)
    final_net.to_csv(f'{sample}.mst.all.csv',sep='\t',header=True,index=False)
    final_net=final_net.groupby(['from','to']).agg('sum').reset_index()
    final_net.to_csv(f'{sample}.mst.merge.csv',sep='\t',header=True,index=False)
    # plot dot
    print_dot(drawpoints,final_net,None)

###########################################################
# entry, usage and parameter
###########################################################
if __name__ == '__main__':
    # usage
    if len(sys.argv) < 3:
        print('CellColocation.py <infile> <sample> [loopnum] [grid_binsize]',flush=True)
        sys.exit(0)
    # default parameters
    mincell = 100 #ignore one celltype if total cell nummber less than 100 
    bootnum = 100
    sample_fac = 0.8 
    binsize=50
    # parse inputs
    label_xyz_file=sys.argv[1]
    sample=sys.argv[2]
    if len(sys.argv)>3:
        bootnum = int(sys.argv[3])
    if len(sys.argv)>4:
        binsize = int(sys.argv[4])
    # call main
    mainpipe(label_xyz_file,sample,bootnum,binsize)
