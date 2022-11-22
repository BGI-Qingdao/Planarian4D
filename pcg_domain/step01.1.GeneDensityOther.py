#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

########################################
# conf
y_lim = 100
min_cell = 100
basic_step = 10 # a little greater than the cell diameter 

########################################
# create the cell count table
class SampleInfo:
    def __init__(self,sample_name):
        self.name = sample_name
        self.pos_file = f'cell_pos_db/{sample_name}.txt'

    def LoadBody(self, y_lim, basic_step):
        #############################
        # load cells
        body = np.loadtxt(self.pos_file)
        bd = pd.DataFrame(columns=['x','y','z']);
        bd['x'] = body[:,0]
        bd['y'] = body[:,1]
        bd['z'] = body[:,2]
        bd['x'] = bd['x'] / basic_step
        bd['x'] = bd['x'].astype(int)
        #############################
        # drop invalid y
        bd = bd[ ( ( bd['y'] < -y_lim) | ( bd['y'] > y_lim ) )].copy()
        bd['z'] = np.ones(len(bd),dtype=int)
        #############################
        # count cell
        bd = bd.groupby(['x']).agg(c=('z','sum')).reset_index()
        bd.to_csv(f'./basic_cache/{self.name}.body_cell_count.other.csv')
        #############################
        # save data
        self.Body = bd
        self.Xmin = int(np.min(bd['x']))
        self.Xmax = int(np.max(bd['x']))
        self.AP_len = int(self.Xmax - self.Xmin + 1)

    def CellCount(self, xmin,xmax):
        return np.sum(self.Body[ ( self.Body['x']>=xmin) & ( self.Body['x']<=xmax) ]['c'])

    def Valid(self, x):
        #return True
        return np.sum(self.Body[self.Body['x']==x]['c']) >= min_cell

########################################
# create the expression table
class GeneInfo:
    def __init__(self,sample_name,gene_name):
        self.name = sample_name
        self.gene = gene_name
        self.exp_file = f'gene_atlas_smes_all/{sample_name}/{gene_name}.txt'

    def LoadExp(self,ylim,basic_step):
        #############################
        # load exp
        cache = np.loadtxt(self.exp_file)
        exp_data = pd.DataFrame(cache);
        exp_data.columns = ['x','y','z','e']
        exp_data['x'] = exp_data['x']/basic_step
        exp_data['x'] = exp_data['x'].astype(int)
        #############################
        # drop invalid y
        exp_data = exp_data[ (exp_data['y'] < -y_lim) | ( exp_data['y'] > y_lim )].copy()
        #############################
        # count exp
        exp_data = exp_data.groupby(['x']).agg({ 'e' :['count','sum']}).reset_index()
        exp_data.columns = ['x','cell_count','exp_sum']       
        exp_data.to_csv(f'./basic_cache/{self.name}_{self.gene}.basic.other.csv')
        #############################
        # save data
        self.exp_data = exp_data
        self.total_cell =int(np.sum(exp_data['cell_count']))
        density_unit = int(self.total_cell/100)
        if density_unit < 10:
            self.density_unit = 10
        else:
            self.density_unit = density_unit
    
    def resetExp(self, xmin, AP_len):
        #print(AP_len,flush=True)
        exp_array = np.zeros((AP_len,2),dtype = float)
        for idx,row in self.exp_data.iterrows():
            idx = int(row['x'] - xmin)
            exp_array[idx,0] = row['cell_count']
            exp_array[idx,1] = row['exp_sum']
        self.xmin = int(xmin)
        self.xmax = int(xmin + AP_len -1)
        self.AP_len = AP_len
        self.exp_array = exp_array
             
    def GetExpIn(self, x_center):
        idx = x_center - self.xmin
        collect_cell = self.exp_array[idx,0] 
        #print(collect_cell,flush=True)
        collect_value = self.exp_array[idx,1] 
        t_min = x_center
        t_max = x_center
        for i in np.arange(1,self.AP_len):
            #print(collect_cell,flush=True)
            if collect_cell > self.density_unit :
                break
            if x_center-i >= self.xmin:
               idx = x_center - self.xmin - i
               collect_cell = collect_cell + self.exp_array[idx,0]
               collect_value = collect_value + self.exp_array[idx,1]
               t_min = x_center-i
            if x_center+i <= self.xmax:
               idx = x_center - self.xmin + i
               collect_cell = collect_cell + self.exp_array[idx,0]
               collect_value = collect_value + self.exp_array[idx,1]
               t_max = x_center+i

        return t_min, t_max, collect_cell,collect_value

#######################################################################################
# main codes
samples = []
APs = []
sample_info = {}
AP_infos = {}
sample_names =  ['WT','0hpa1','0hpa2','12hpa1','12hpa2','36hpa1','36hpa2','3dpa1','3dpa2',
                '5dpa1','5dpa2','7dpa1','7dpa2','10dpa1','10dpa2','14dpa1','14dpa2' ] 
pcgs = np.loadtxt('target_genes',dtype=str)

for x in sample_names:
    si = SampleInfo(x)
    si.LoadBody(y_lim, basic_step)
    sample = [x] * si.AP_len
    AP = np.arange(si.Xmin,si.Xmax+1,1)
    samples.append(np.array(sample))
    AP_infos[x] = AP
    APs.append(AP) 
    sample_info[x] = si

raw_data = pd.DataFrame(columns=['sample','AP','t_min','t_max','exp_cell_count','exp_sum','cell_count','density1k'])
raw_data['sample'] = np.concatenate(samples,axis=0)
raw_data['AP'] = np.concatenate(APs,axis=0)

all_gene_data =  pd.DataFrame(columns=['sample','AP','t_min','t_max','exp_cell_count','exp_sum','cell_count','density1k','gene'])
for smes in pcgs:
    t_min_l = []
    t_max_l = []
    cc_l = []
    cv_l = []
    ed_l = []
    current_gene = raw_data.copy()
    for x in sample_names:
        si = sample_info[x]
        ei = GeneInfo(x,smes)
        ei.LoadExp(y_lim, basic_step)
        ei.resetExp(si.Xmin, si.AP_len)
        cache = np.zeros((len(AP_infos[x]),6))
        for i in range(len(AP_infos[x])):
            x_center = AP_infos[x][i]
            if si.Valid(x_center) == False:
                continue
            ti, tm, cc, ev = ei.GetExpIn(x_center)
            cache[i,0] = ti
            cache[i,1] = tm
            cache[i,2] = cc
            cache[i,3] = ev
            cell_count = si.CellCount(ti,tm)
            cache[i,4] = cell_count
            cache[i,5] = float(ev * 1000 )/float(cell_count)
        current_gene.loc[(current_gene['sample']== x) ,'t_min'] = cache[:,0]
        current_gene.loc[(current_gene['sample']== x) ,'t_max'] = cache[:,1]
        current_gene.loc[(current_gene['sample']== x) ,'exp_cell_count'] = cache[:,2]
        current_gene.loc[(current_gene['sample']== x) ,'exp_sum'] = cache[:,3]
        current_gene.loc[(current_gene['sample']== x) ,'cell_count'] = cache[:,4]
        current_gene.loc[(current_gene['sample']== x) ,'density1k'] = cache[:,5]
    current_gene['gene'] = [smes] * len(current_gene)
    current_gene = current_gene.drop(index = current_gene[ current_gene['cell_count'] == 0 ].index[:])
    if len(all_gene_data) == 0 :
        all_gene_data = current_gene
    else:
        all_gene_data = pd.concat([all_gene_data,current_gene],ignore_index=True)

all_gene_data.to_csv('raw_all.ml100.other.csv')
