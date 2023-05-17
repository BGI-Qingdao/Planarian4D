#!/usr/bin/env python3

###########################################################
# imports
###########################################################
import sys
import json
import numpy as np
import pandas as pd

###########################################################
# functions
###########################################################

def get_trackEM_forward_affine(param :str) -> np.matrix:
    """
        handle '-0.010963829,-0.999939895,0.999939895,-0.010963829,-129.2603788,1664.628308'

        @return affine matrix.
    """
    affine = np.zeros((3,3))
    in_data = np.array(param.split(',')).astype(float).reshape((3,2))
    affine[0:2,:]=in_data.T
    affine[2] = [0,0,1]
    return np.matrix(affine)

def mask_gemc(gemc,mask):
    ####
    # for 5dap2-slice12 
    ymax = np.max(gemc['y'])+1
    xmax = np.max(gemc['x'])+1
    yl,xl = mask.shape
    if ymax>yl or xmax>xl:
        print(f'{ymax},{xmax},{yl},{xl}')
        print('Waring ... heatmap not fit mask!!!!!!!!!!!!!!')
        new_mask = np.zeros((max(ymax,yl),max(xmax,xl)),dtype=int)
        y,x = np.nonzero(mask)
        new_mask[y,x]=1
        mask=new_mask
    gemc['mask'] = mask[gemc['y'],gemc['x']]
    gemc = gemc[gemc['mask'] != 0].copy()
    gemc.drop(columns=['mask'],inplace=True)
    return gemc
 
###########################################################
# main class
###########################################################
class Cache:
    def __init__(self,prefix):
        self.prefix = prefix
        self.tags = []
        self.masks = {}
        self.sids = {}
        self.affines = {}
        self.minx = {}
        self.miny = {}

    def load_conf(self):
        prefix = self.prefix
        conf = json.load(open(f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/s05.3D_affine/s02.3d_affine_conf/{prefix}.json'))
        for clist in conf:
            slice_id = clist[0]
            self.tags.append(slice_id)
            self.sids[slice_id] = int(clist[1])
            self.masks[slice_id] = clist[2]
            self.affines[slice_id] = get_trackEM_forward_affine(clist[3])

    def load_minxy(self):
        conf = json.load(open('/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/s05.3D_affine/s02.3d_affine_conf/GEMC_ROI.txt'))
        for clist in conf:
            slice_id = clist[0]
            indv = clist[1]
            if indv == self.prefix:
                self.minx[slice_id] = clist[2]
                self.miny[slice_id] = clist[3]

    
    def load_gemc(self, sid):
        prefix = self.prefix
        gemc_file =  f'/dellfsqd2/ST_OCEAN/USER/guolidong/wochong_rna_202201/s04.final_gemc/{prefix}/{prefix}_{sid}.gemc'
        #geneID  x       y       MIDCounts       cell
        data = pd.read_csv(gemc_file,sep='\t',header=0)
        ll = len(data)
        minx = self.minx[sid]
        miny = self.miny[sid]
        print(f'load #RNA {ll} from {gemc_file}, with x>{minx}, y>{miny}',flush=True)
        data['x'] = data['x'] - minx
        data['y'] = data['y'] - miny
        return data
    
    def load_mask(self,tag):
        mask_file = self.masks[tag]
        data = np.loadtxt(mask_file,delimiter=' ',dtype=int)
        return data

    def affine(self, gemc, tag):
        gemc['3d_z'] = np.ones(len(gemc))
        affine_matrix = self.affines[tag]
        raw_coord = gemc[['x','y','3d_z']].to_numpy()
        new_coord = np.dot(affine_matrix,raw_coord.T)
        new_coord = new_coord.T
        gemc['3d_x']= new_coord[:,0]
        gemc['3d_y']= new_coord[:,1]
        gemc['3d_z']= new_coord[:,2]
        return gemc

    def set_xyz_3dxyz(self,data,sid):
        minx = self.minx[sid]
        miny = self.miny[sid]
        data['x'] = data['x'] + minx
        data['y'] = data['y'] + miny
        data['3d_z'] = data['3d_z']*sid
        data['3d_x'] = data['3d_x']*10
        data['3d_y'] = data['3d_y']*10
        data['3d_z'] = data['3d_z']*140
        data['3d_x'] = data['3d_x'].astype(int)
        data['3d_y'] = data['3d_y'].astype(int)
        data['3d_z'] = data['3d_z'].astype(int)
        return data       

    def affineMain(self):
        self.load_conf()
        self.load_minxy()
        for tag in self.tags:
            sid = self.sids[tag]
            gemc = self.load_gemc(sid)
            mask = self.load_mask(tag)
            gemc =      mask_gemc(gemc,mask)
            gemc = self.affine(gemc,tag)
            gemc = self.set_xyz_3dxyz(gemc,sid)
            gemc.to_csv(f'{self.prefix}.{tag}.gemc_3d',sep='\t',header=True,index=False)

if __name__ == '__main__':
    prefix = sys.argv[1]
    cache = Cache(prefix)
    cache.affineMain()
