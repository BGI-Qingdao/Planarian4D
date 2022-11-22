#!/usr/bin/env python3
import sys
import getopt
import numpy as np
import pandas as pd
from skimage import io as skio

#################################################
# BinConf
# 
class BinConf:
    def __init__(self,view,binsize):
        self.gene_binsize = binsize
        if view == 'APML' :
            minBorderBinsize = 10
        else :
            minBorderBinsize = 20 # DV may contain 20 micron gaps
        if binsize < minBorderBinsize :
            self.body_binsize = minBorderBinsize
            self.body_scale = int(minBorderBinsize/binsize) #please use binsize = 5/10/20 or greater value !
        else :
            self.body_binsize = binsize
            self.body_scale = 1

    def geneBinsize(self):
        return self.gene_binsize

    def bodyBinConf(self):
        return self.body_binsize , self.body_scale

#################################################
# BorderDetect 
# 
class BorderDetect:
    def __init__(self,x,y):
        self.x = x 
        self.y = y

    def Border(self):
        xmin = np.min(self.x)
        xmax = np.max(self.x)
        ymin = np.min(self.y)
        ymax = np.max(self.y)
        self.x = self.x -xmin + 5  #left 5 pixel of margin
        self.y = self.y -ymin + 5  #left 5 pixel of margin
        height = int(ymax-ymin+10) #right 5 pixel of margin
        width = int(xmax-xmin+10)  #right 5 pixel of margin
        mask = np.zeros((height,width),dtype=int)
        mask[self.y,self.x] = 1
        # open and close to remove noise and fill holes
        from scipy import ndimage
        mask = ndimage.binary_closing(mask).astype(int)
        mask = ndimage.binary_opening(mask).astype(int)
        
        for y in range(0,height):
            for x in range(1,width-1):
                if mask[y,x-1]==0 and mask[y,x]==0 and mask[y,x+1]==1:
                    mask[y,x] = 255
            for x in range(width-2,1,-1):
                if mask[y,x+1]==0 and mask[y,x]==0 and mask[y,x-1]==1:
                    mask[y,x] = 255
        for x in range(0,width):
            for y in range(1,height-1):
                if mask[y-1,x]==0 and mask[y,x]==0 and mask[y+1,x]==1:
                    mask[y,x] = 255
            for y in range(height-2,1,-1):
                if mask[y+1,x]==0 and mask[y,x]==0 and mask[y-1,x]==1:
                    mask[y,x] = 255
        mask[self.y,self.x] = 0
        (y_idx,x_idx) = np.nonzero(mask)
        y_idx = y_idx + ymin - 5 #switch back to raw coord
        x_idx = x_idx + xmin - 5 #switch back to raw coord
        #xy = np.zeros((len(y_idx),2))
        #xy[:,0] = x_idx
        #xy[:,1] = y_idx
        #import alphashape
        #hull = alphashape.alphashape(xy,alpha=0.001)
        #hull_xy = np.array(list(zip(*hull.exterior.xy)))
        #hull_xy = hull_xy.astype(int)
        #return hull_xy[:,1] , hull_xy[:,0]
        return y_idx,x_idx

#################################################
# ROIManager
#
class ROIManager:
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
    
    ###################################
    # always filter nax before min
    def filterDFMax(self,df):
        if self.zmax != None:
            df = df[df['z']<=self.zmax].copy()
        if self.xmax != None:
            df = df[df['x']<=self.xmax].copy()
        if self.ymax != None:
            df = df[df['y']<=self.ymax].copy()
        return df

    ###################################
    # always filter nax before min
    def filterAndResetDFMin(self,df):
        if self.xmin == None:
            self.xmin = np.min(df['x'])
        if self.ymin == None:
            self.ymin = np.min(df['y'])
        if self.zmin == None:
            self.zmin = np.min(df['z'])
        df['x'] = df['x'] - self.xmin
        df = df[df['x']>=0].copy()
        df['y'] = df['y'] - self.ymin
        df = df[df['y']>=0].copy()
        df['z'] = df['z'] - self.zmin
        df = df[df['z']>=0].copy()
        return df
            

#################################################
# BodyInfo 
# 
class BodyInfo:
    def __init__(self,bin_border,bin_draw_scale):
        self.bin_border = bin_border
        self.bin_draw_scale = bin_draw_scale

    def loadAllPoints(self,xyz_txt, roi):
        body = np.loadtxt(xyz_txt)
        bd = pd.DataFrame(columns=['x','y','z']);
        #format x axis by roi first
        bd['x'] = body[:,0]
        bd['z']= body[:,2]
        bd['y']= body[:,1]
        bd = roi.filterDFMax(bd)
        bd = roi.filterAndResetDFMin(bd)
        #bin x axis now
        bd['x']= bd['x']/self.bin_border
        bd['x']= bd['x'].astype(int)
        bd['x']=bd['x']+1
        #bin y axis now
        bd['y']= bd['y']/self.bin_border
        bd['y']= bd['y'].astype(int)
        bd['y']= bd['y']+1
        #bin z axis now
        bd['z']= bd['z']/self.bin_border
        bd['z']= bd['z'].astype(int)
        bd['z']= bd['z']+1
        # save final body
        self.body = bd

    #################################
    #  AP = x axis , ML = y axis
    #
    def calcAPML_border(self):              
        # agg all z coordinate
        bd = self.body.groupby(['x', 'y']).agg(z=('z','max')).reset_index()
        height = int(np.max(bd['y'])+2) #right 1 pixel margin
        width = int(np.max(bd['x'])+2)  #right 1 pixel margin
        # get basic infos
        self.APML_W  = width * self.bin_draw_scale
        self.APML_H  = height * self.bin_draw_scale
        self.APML_points_num = len(bd)*self.bin_draw_scale*self.bin_draw_scale
        # get border dash points
        graph_matrix = np.zeros((height,width),dtype='uint8')
        graph_matrix[bd['y'],bd['x']]=1
        ( body_y, body_x ) = np.nonzero(graph_matrix)
        y_idx,x_idx = BorderDetect(body_x,body_y).Border()
        # save final border
        self.APML_x_idx = x_idx*self.bin_draw_scale
        self.APML_y_idx = y_idx*self.bin_draw_scale

    def getAPML_num_points(self):
        return self.APML_points_num

    def getAPML_WH(self):
        return self.APML_W,self.APML_H

    def getAPML_border(self):   
        return self.APML_x_idx , self.APML_y_idx           
    
    #################################
    #  AP = x axis , DV = y axis
    #
    def calcAPDV_border(self):              
        # agg all y coordinate
        bd = self.body.groupby(['x', 'z']).agg(y=('y','max')).reset_index()
        height = int(np.max(bd['z'])+2)
        width = int(np.max(bd['x'])+2)
        # get basic infos
        self.APDV_W  = width * self.bin_draw_scale
        self.APDV_H  = height * self.bin_draw_scale
        self.APDV_points_num = len(bd)*self.bin_draw_scale*self.bin_draw_scale
        # get border dash points
        graph_matrix = np.zeros((height,width),dtype='uint8')
        graph_matrix[bd['z'],bd['x']]=1
        ( body_y, body_x ) = np.nonzero(graph_matrix)
        y_idx,x_idx = BorderDetect(body_x,body_y).Border()
        # save final border
        self.APDV_x_idx = x_idx*self.bin_draw_scale
        self.APDV_y_idx = y_idx*self.bin_draw_scale

    def getAPDV_num_points(self):
        return self.APDV_points_num

    def getAPDV_WH(self):
        return self.APDV_W, self.APDV_H

    def getAPDV_border(self):
        return self.APDV_x_idx , self.APDV_y_idx

    #################################
    #  ML = x axis , DV = x axis
    #
    def calcMLDV_border(self):              
        # agg all y coordinate
        bd = self.body.groupby(['y', 'z']).agg(x=('x','max')).reset_index()
        height = int(np.max(bd['y'])+2)
        width = int(np.max(bd['z'])+2)
        # get basic infos
        self.MLDV_W  = width * self.bin_draw_scale
        self.MLDV_H  = height * self.bin_draw_scale
        self.MLDV_points_num = len(bd)*self.bin_draw_scale*self.bin_draw_scale
        # get border dash points
        graph_matrix = np.zeros((height,width),dtype='uint8')
        graph_matrix[bd['y'],bd['z']]=1
        ( body_y, body_x ) = np.nonzero(graph_matrix)
        y_idx,x_idx = BorderDetect(body_x,body_y).Border()
        # save final border
        self.MLDV_x_idx = x_idx*self.bin_draw_scale
        self.MLDV_y_idx = y_idx*self.bin_draw_scale

    def getMLDV_num_points(self):
        return self.MLDV_points_num

    def getMLDV_WH(self):
        return self.MLDV_W, self.MLDV_H

    def getMLDV_border(self):   
        return self.MLDV_x_idx , self.MLDV_y_idx           

class Gene3D:
    def __init__(self,binsize):
        self.binsize = binsize

    def loadExpr(self,gene_txt_list,roi):
        a = ''
        for gene_txt in gene_txt_list :  #hanle multi files here
            aa = np.loadtxt(gene_txt)
            if len(a) == 0:
               a = aa
            else :
               a = np.concatenate((a, aa), axis=0)
            
        if len(a) == 0 :
           self.valid = False   
           return
        show_data = pd.DataFrame(columns=['x','y','z','value'])
        show_data['x'] = a[:,0]
        show_data['y'] = a[:,1]
        show_data['z'] = a[:,2]
        show_data['value'] = a[:,3]
        show_data = roi.filterDFMax(show_data)
        show_data = roi.filterAndResetDFMin(show_data)
        if len(show_data) <1 :
            self.valid = False
            return
        self.valid = True
        self.gene_expr = show_data

    def getMIR_APML(self):
        show_data = self.gene_expr.copy()
        show_data['x'] = show_data['x']/self.binsize
        show_data['x'] = show_data['x'].astype(int)
        show_data['y'] = show_data['y']/self.binsize
        show_data['y'] = show_data['y'].astype(int)
        show_data = show_data.groupby(['x', 'y']).agg(value=('value', 'max')).reset_index()
        return show_data

    def getMIR_APDV(self):
        show_data = self.gene_expr.copy()
        show_data['x'] = show_data['x']/self.binsize
        show_data['x'] = show_data['x'].astype(int)
        show_data['z'] = show_data['z']/self.binsize
        show_data['z'] = show_data['z'].astype(int)
        show_data = show_data.groupby(['x', 'z']).agg(value=('value', 'max')).reset_index()
        return show_data

    def getMIR_MLDV(self):
        show_data = self.gene_expr.copy()
        show_data['y'] = show_data['y']/self.binsize
        show_data['y'] = show_data['y'].astype(int)
        show_data['z'] = show_data['z']/self.binsize
        show_data['z'] = show_data['z'].astype(int)
        show_data = show_data.groupby(['y', 'z']).agg(value=('value', 'max')).reset_index()
        return show_data
   
def GetBodyInfo(body_txt,binconf,roi):
    body_binsize, body_scale = binconf.bodyBinConf()
    body_info = BodyInfo(body_binsize,body_scale)
    body_info.loadAllPoints(body_txt,roi)
    return body_info

def GetBackground(view,body_info):
    if view == "APML" :
        body_info.calcAPML_border()
        W,H = body_info.getAPML_WH()
        draw_array = np.zeros((H,W,3),dtype='int')
        xids,yids = body_info.getAPML_border()
        # draw background
        draw_array[yids,xids,:] = 255
        # draw scale bar
        #draw_array[[H-10,H-9,H-8],W-15:W-5,:]=255
        #draw_array[H-20:H-10,W-15:W-12,:]=255
        draw_array[H-50:H-10,W-16:W-10,:]=255
        return draw_array
    elif view == "APDV" :
        body_info.calcAPDV_border()
        W,H = body_info.getAPDV_WH()
        draw_array = np.zeros((H,W,3),dtype='int')
        xids,yids = body_info.getAPDV_border()
        # draw background
        draw_array[yids,xids,:] = 255
        # draw scale bar
        #draw_array[[H-10,H-9,H-8],W-15:W-5,:]=255
        draw_array[H-20:H-10,W-15:W-12,:]=255
        return draw_array
        #return None
    elif view == "MLDV" :
        body_info.calcMLDV_border()
        W,H = body_info.getMLDV_WH()
        draw_array = np.zeros((H,W,3),dtype='int')
        xids,yids = body_info.getMLDV_border()
        # draw background
        draw_array[yids,xids,:] = 255
        # draw scale bar
        #draw_array[[H-10,H-9,H-8],W-15:W-5,:]=255
        draw_array[H-20:H-10,W-15:W-12,:]=255
        return draw_array

def GetGeneExpr(gene_txt,binconf,roi):
    gene_expr = Gene3D(binconf.geneBinsize())
    gene_expr.loadExpr(gene_txt,roi)
    return gene_expr

def FISH_scale(num_points, panel_expr):
    from sklearn.preprocessing import QuantileTransformer
    min_value = np.min(panel_expr)
    max_value = np.max(panel_expr)
    #print(f'min = {min_value}; max={max_value}',flush=True)
    num_light_panel=len(panel_expr)
    #temp_list = np.zeros(num_points)
    #temp_list[0:num_light_panel]=panel_expr
    #qt = QuantileTransformer(output_distribution='uniform', random_state=0)
    #temp_list_scale = qt.fit_transform(temp_list.reshape(-1, 1))
    #temp_list_scale = temp_list_scale.reshape(-1)
    #ret_data = temp_list_scale[0:num_light_panel]
    ret_data = panel_expr / max_value
    ret_data = ret_data *255
    #ret_data = np.power((np.ones(len(ret_data))*2),ret_data);
    ret_data [ ret_data >255 ] = 255
    ret_data = ret_data.astype(int)
    return ret_data
    
def DrawSingleFISH_APML( body_info, expr, channel_id):
    W,H = body_info.getAPML_WH()
    draw_array = np.zeros((H,W,3),dtype='uint8')
    APML_expr = expr.getMIR_APML()
    APML_expr['y'] = APML_expr['y']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    APML_expr['x'] = APML_expr['x']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    draw_expr =  FISH_scale(body_info.getAPML_num_points(),APML_expr['value'])
    if channel_id == 0: 
        draw_array[APML_expr['y'],APML_expr['x'],0] = draw_expr
        draw_array[APML_expr['y'],APML_expr['x'],2] = draw_expr
    elif channel_id == 2: 
        draw_array[APML_expr['y'],APML_expr['x'],0] = draw_expr
        draw_array[APML_expr['y'],APML_expr['x'],1] = draw_expr
    else:
        draw_array[APML_expr['y'],APML_expr['x'],channel_id] = draw_expr
    return draw_array 

def DrawSingleFISH_APDV(body_info, expr, channel_id):
    W,H = body_info.getAPDV_WH()
    draw_array = np.zeros((H,W,3),dtype='uint8')
    APDV_expr = expr.getMIR_APDV()
    APDV_expr['x'] = APDV_expr['x']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    APDV_expr['z'] = APDV_expr['z']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    draw_expr = FISH_scale(body_info.getAPDV_num_points(),APDV_expr['value'])
    if channel_id == 0: 
        draw_array[APDV_expr['z'],APDV_expr['x'],0] = draw_expr
        draw_array[APDV_expr['z'],APDV_expr['x'],2] = draw_expr
    elif channel_id == 2: 
        draw_array[APDV_expr['z'],APDV_expr['x'],0] = draw_expr
        draw_array[APDV_expr['z'],APDV_expr['x'],1] = draw_expr
    else:
        draw_array[APDV_expr['z'],APDV_expr['x'],channel_id] = draw_expr
    return draw_array 

def DrawSingleFISH_DVML(body_info, expr, channel_id):
    W,H = body_info.getMLDV_WH()
    draw_array = np.zeros((H,W,3),dtype='uint8')
    MLDV_expr = expr.getMIR_MLDV()
    MLDV_expr['y'] = MLDV_expr['y']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    MLDV_expr['z'] = MLDV_expr['z']+body_info.bin_draw_scale #shift the 1 pixel margin for border
    draw_expr = FISH_scale(body_info.getMLDV_num_points(),MLDV_expr['value'])
    if channel_id == 0: 
        draw_array[MLDV_expr['y'],MLDV_expr['z'],0] = draw_expr
        draw_array[MLDV_expr['y'],MLDV_expr['z'],2] = draw_expr
    elif channel_id == 2: 
        draw_array[MLDV_expr['y'],MLDV_expr['z'],0] = draw_expr
        draw_array[MLDV_expr['y'],MLDV_expr['z'],1] = draw_expr
    else:
        draw_array[MLDV_expr['y'],MLDV_expr['z'],channel_id] = draw_expr
    return draw_array 

def DrawSingleFISH(view, body_info, gene_expr, channel_id):
    if view == "APML" :
        return DrawSingleFISH_APML(body_info, gene_expr, channel_id)
    elif view == "APDV":
        return DrawSingleFISH_APDV(body_info, gene_expr, channel_id)
    elif view == "MLDV":
        return DrawSingleFISH_DVML(body_info, gene_expr, channel_id)

############################################################################
# Usage
#
def FISH_like_usage():
    print("""
Usage : FISH_3c.py -i <individual.txt> \\
                   -o <output prefix>  \\
                   -g [gene.txt that draw in Green(#00ff00) channel] \\
                   -m [gene.txt that draw in Magenta(#ff00ff) channel]   \\
                   -y [gene.txt that draw in Yellow(#ffff00) channel]  \\
                   --view [default APML, must be APML/APDV/MLDV] \\
                   --xmin [default None] \\
                   --ymin [default None] \\
                   --zmin [default None] \\
                   --xmax [default None] \\
                   --ymax [default None] \\
                   --zmax [default None] \\
                   --binsize [default 5] 
Example :
    python3 FISH_3c.py -i WT.txt -o WT_wnt1 -g wnt1.txt 
Notice :
    please at least assign one of -g, -m, or -y.
    -g/-m/-y can be assigned many times.
""")


def saveimage(fname ,draw_matrix):
    used_draw_matrix = draw_matrix.copy()
    used_draw_matrix[draw_matrix>255] = 255
    new_draw_matrix = np.zeros(used_draw_matrix.shape , dtype='uint8') 
    #used_draw_matrix = used_draw_matrix.astype('uint8')
    new_draw_matrix[:,:,:] = used_draw_matrix[:,:,:];
    skio.imsave(fname,new_draw_matrix)
    

def FISH_like_main(argv:[]):
    ###############################################################################
    # Default values
    indv = ''
    prefix = ''
    g_gene = []
    m_gene = []
    y_gene = []
    view = 'APML'
    xmin = ymin = zmin = None
    xmax = ymax = zmax = None
    binsize = 5

    ###############################################################################
    # Parse the arguments
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:g:y:",["help","view=","xmin=","ymin=","zmin=","xmax=","ymax=","zmax=","binsize="])
    except getopt.GetoptError:
        FISH_like_usage()
        sys.exit(2)
    for opt, arg in opts:
        print(f' {opt} {arg} ')
        if opt in ('-h' ,'--help'):
            FISH_like_usage()
            sys.exit(0)
        elif opt in ("-i"):
            indv = arg
        elif opt in ("-o"):
            prefix = arg
        elif opt == "-m" :
            m_gene.append(arg)
        elif opt == "-y" :
            y_gene.append(arg)
        elif opt == "-g" :
            g_gene.append(arg)
        elif opt == "--xmin":
            xmin = int(arg)
        elif opt == "--ymin":
            ymin = int(arg)
        elif opt == "--zmin":
            zmin = int(arg)
        elif opt == "--xmax":
            xmax = int(arg)
        elif opt == "--ymax":
            ymax = int(arg)
        elif opt == "--zmax":
            zmax = int(arg)
        elif opt == "--view":
            view = arg
        elif opt == "--binsize":
            binsize = int(arg)
    
    ###############################################################################
    # Sanity check
    if indv == "" or ( len(g_gene) == 0 and len(m_gene) == 0 and len(y_gene) == 0 ) or prefix == "":
        FISH_like_usage()
        sys.exit(3)
    roi = ROIManager(xmin,xmax,ymin,ymax,zmin,zmax)
    binconf = BinConf(view,binsize)
    print(f"the drawing view : {view}")
    ###############################################################################
    # Load the body points 
    print('Loading body now ...',flush=True)
    body_info = GetBodyInfo(indv,binconf,roi)
    ###############################################################################
    # Load the gene expr points and draw
    print('Loading gene expression now ...',flush=True)
    m_gene_expr = g_gene_expr = y_gene_expr = None 
    channel = 0
    if len(g_gene) != 0 :
        g_gene_expr = GetGeneExpr(g_gene,binconf,roi) 
        if g_gene_expr.valid :
            channel = channel + 1
        else :
            g_gene_expr = None
    if m_gene != "" :
        m_gene_expr = GetGeneExpr(m_gene,binconf,roi) 
        if m_gene_expr.valid :
            channel = channel + 1
        else :
            m_gene_expr = None
    if len(y_gene) != 0  :
        y_gene_expr = GetGeneExpr(y_gene,binconf,roi) 
        if y_gene_expr.valid :
            channel = channel + 1
        else :
            y_gene_expr = None

    ###############################################################################
    # get sample border
    background = GetBackground(view,body_info)
    ###############################################################################
    # Draw single FISH
    r_image = b_image = g_image = None
    if m_gene_expr != None :
        r_image = DrawSingleFISH(view,body_info,m_gene_expr,0)
        if channel == 1 :
            print('Draw magenta-channel-FISH now ...',flush=True)
            draw_image = background + r_image
            saveimage(f'{prefix}.magenta.MIR.png',draw_image)
            #skio.imsave(f'{prefix}.magenta.MIR.tif',draw_image.astype('uint8'))
    if g_gene_expr != None :
        g_image = DrawSingleFISH(view,body_info,g_gene_expr,1)
        if channel == 1 :
            print('Draw green-channel-FISH now ...',flush=True)
            draw_image = background + g_image 
            saveimage(f'{prefix}.green.MIR.png',draw_image)
            #skio.imsave(f'{prefix}.green.MIR.tif',draw_image.astype('uint8'))
    if y_gene_expr != None :
        b_image = DrawSingleFISH(view,body_info,y_gene_expr,2)
        if channel == 1 :
            print('Draw yellow-channel-FISH now ...',flush=True)
            draw_image = background + b_image 
            saveimage(f'{prefix}.yellow.MIR.png',draw_image)
            #skio.imsave(f'{prefix}.yellow.MIR.tif',draw_image.astype('uint8'))

    ###############################################################################
    # Draw multi-FISH
    if channel > 1:
        print('Draw multi-FISH now ...',flush=True)
        if m_gene_expr != None :
            draw_image = background + r_image
        if g_gene_expr != None :
            draw_image = draw_image + g_image
        if y_gene_expr != None :
            draw_image = draw_image + b_image
        saveimage(f'{prefix}.multi.MIR.png',draw_image)
        #skio.imsave(f'{prefix}.multi.MIR.tif',draw_image)
    
    ###############################################################################
    # Done
    print('__ALL DONE__',flush=True)

if __name__ == "__main__":
    FISH_like_main(sys.argv[1:])
