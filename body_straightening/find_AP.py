#!/usr/bin/env python3
import pandas as pd
import numpy as np
from skimage import io
import scipy.ndimage as nd
import scipy.linalg as sl



def polyfix(x, y, n, xfix, yfix, xder='', dydx=''):
    """ Fits polynomial p with degree n to data,
        but specify value at specific points
    """

    nfit = len(x)
    if len(y) != nfit:
        raise ValueError('x and y must have the same size')
    nfix = len(xfix)
    if len(yfix) != nfix:
        raise ValueError('xfit adn yfit must have the same size')
    x = np.vstack(x)
    y = np.vstack(y)
    xfix = np.vstack(xfix)
    yfix = np.vstack(yfix)

    # if derivatives are specified:
    if xder != '':
        try:
            len(xder)
            xder = xder
            dydx = dydx
        except:
            xder = np.array([xder])
            dydx = np.array([dydx])

    else:
        xder = []
        dydx = []

    nder = len(xder)
    if len(dydx) != nder:
        raise ValueError('xder and dydx must have same size')

    nspec = nfix + nder

    ##specval = np.vstack((yfix, dydx))
    specval = yfix
    # first find A and pc such that A*pc = specval
    A = np.zeros((nspec, n+1))
    # specified y values
    for i in range(n+1):
        A[:nfix, i] = np.hstack(np.ones((nfix, 1)) * xfix**(n+1-(i+1)))
    if nder > 0:
        for i in range(n):
            A[nfix:nder+nfix, i] = ((n-(i+1)+1) * np.ones((nder, 1)) * xder**(n-(i+1))).flatten()
    if nfix > 0:
        lastcol = n+1
        nmin = nspec - 1
    else:
        lastcol = n
        nmin = nspec

    if n < nmin:
        raise ValueError('Polynomial degree too low, cannot match all constraints')
    # find unique polynomial of degree nmin that fits the constraints
    firstcol = n-nmin#+1
    pc0 = np.linalg.solve(A[:, firstcol:lastcol], specval)
    pc = np.zeros((n+1, 1))
    pc[firstcol:lastcol] = pc0

    X = np.zeros((nfit, n+1))
    for i in range(n+1):
        X[:, i] = (np.ones((nfit, 1)) * x**(n+1-(i+1))).flatten()

    yfit = y - np.polyval(pc, x)

    B = sl.null_space(A)
    z = np.linalg.lstsq(X @ B, yfit,rcond=None)[0]
    if len(z) == 0:
        z = z[0]
        p0 = B*z
    else:
        p0 = B@z
    p = np.transpose(p0) + np.transpose(pc)
    return p


def poly_fit3(x, y, w, h, x_max, x_min, y_max, y_min, xf, yf, exponential_number=3):
    f1 = polyfix(x.copy(), y.copy(),exponential_number, xf, yf, )
    print(f1[0])
    xvals = np.array(range(0,w))
    yvals = np.polyval(f1[0], xvals)
    line = np.zeros((h,w))
    line_df = pd.DataFrame(columns=['x','y'])
    line_df['x'] = xvals
    line_df['x'] = (line_df['x']+0.5).astype('int')
    line_df['y'] = yvals
    line_df['y'] = (line_df['y']+0.5).astype('int')
    line_df = line_df[line_df['y']<y_max]
    line_df = line_df[line_df['x']<x_max]
    line_df = line_df[line_df['y']>y_min]
    line_df = line_df[line_df['x']>x_min]
    line[line_df['y'],line_df['x']]=255
    return line,line_df,f1

def poly_fit2(x, y,w, h, x_max, x_min, y_max, y_min, exponential_number=3):
    f1 = np.polyfit(x[1:len(x) - 1], y[1:len(y) - 1], exponential_number)
    xvals = np.array(range(0, w))
    yvals = np.polyval(f1, xvals)
    line = np.zeros((h,w))
    line_df = pd.DataFrame(columns=['x','y'])
    line_df['x'] = xvals
    line_df['x'] = (line_df['x']+0.5).astype('int')
    line_df['y'] = yvals
    line_df['y'] = (line_df['y']+0.5).astype('int')
    line_df = line_df[line_df['y']<y_max]
    line_df = line_df[line_df['x']<x_max]
    line_df = line_df[line_df['y']>y_min]
    line_df = line_df[line_df['x']>x_min]
    line[line_df['y'],line_df['x']]=255
    return line,line_df

def RGB2wh(img):
    new_data = np.zeros((img.shape[0], img.shape[1]), dtype=int)
    new_data = new_data + img[:, :, 0]
    new_data = new_data + img[:, :, 1]
    new_data = new_data + img[:, :, 2]
    new_data = (new_data) / 3
    new_data = new_data.astype('uint8')
    return new_data


def colour_overlap2(img1,img2):
    slice1 = img1.copy()
    slice2 = img2.copy()
    slice1[slice1 != 0] = 1
    slice2[slice1 == 1] = [255, 255, 255]
    return slice2


p = pd.DataFrame(columns = ['id','k0','k1','k2','k3','x_begin','x_end','y_begin','y_end'])

for i in 'WT','0hpa1','0hpa2','3dpa1','3dpa2','5dpa1','5dpa2','7dpa1','7dpa2','10dpa1','10dpa2','12hpa1','12hpa2','14dpa1','14dpa2','36hpa1','36hpa2':

    test = io.imread(f'input/slit/{i}_SMED30023930.green.MIR.tif')
    rawg = test.copy()
    test = RGB2wh(test)
    raw = io.imread(f'input/AMP/{i}.multi.MIR.tif')
    raw = raw[0:raw.shape[0]-10,0:raw.shape[1]-10]
    h = test.shape[0]
    w = test.shape[1]

    #domain
    yt,xt = np.where(test==255)
    yt_max = yt.max()
    yt_min = yt.min()
    xt_max = xt.max()
    xt_min = xt.min()

    test[test==255]=0
    test[test!=0]=255

    y, x = np.where(test != 0)
    #pos = pd.DataFrame()
    #pos['x'] = x
    #pos['y'] = y
    #pos.to_csv(f'output/input_pos/{i}_pos.txt', sep='\t', header=True, index=False)

    constrain = pd.read_csv(f'input/constrain.txt', sep="\t")

    xf = [constrain[constrain['id']==i]['x1'],constrain[constrain['id']==i]['x2']]
    yf = [constrain[constrain['id']==i]['y1'],constrain[constrain['id']==i]['y2']]


    #line1,line_df1 = poly_fit(x,y,w,h,xt_max,xt_min,yt_max,yt_min,xf,yf,3)
    #line2,line_df2 = poly_fit2(x,y,w,h,xt_max,xt_min,yt_max,yt_min,4)
    line,line_df,f1 = poly_fit3(x,y,w,h,xt_max,xt_min,yt_max,yt_min,xf,yf,3)

    #colour1 = colour_overlap(line1,raw)
    #colour2 = colour_overlap2(line2,raw)
    colour = colour_overlap2(line,raw)
    colourg = colour_overlap2(line, rawg)

    #colour1 = colour1.astype('uint8')
    #io.imsave(f'C:/Users/jizhi1/Desktop/correct/result/{i}_colour_c.tif', colour1)
    #colour2 = colour2.astype('uint8')
    #io.imsave(f'output/slit_line/{i}_line_image.tif', colour2)
    #colour = colour.astype('uint8')
    io.imsave(f'output/AMP_line/{i}_line_image.tif', colour)
    io.imsave(f'output/slit_line/{i}_line_image.tif', colourg)

    #line1 = line1.astype('uint8')
    #io.imsave(f'C:/Users/jizhi1/Desktop/correct/result/{i}_line_c.tif', line1)
    #line2 = line2.astype('uint8')
    #io.imsave(f'C:/Users/jizhi1/Desktop/correct/AMP_line/line/{i}_line.tif', line2)
    #line = line.astype('uint8')
    #line_df.to_csv(f'output/line_pos/{i}_line_pos.txt',sep='\t',header=True,index=False)
    p = p.append({'id':i,'k0':f1[0][0],'k1':f1[0][1],'k2':f1[0][2],'k3':f1[0][3],'x_begin':line_df['x'].min(),'x_end':line_df['x'].max(),'y_begin':line_df['y'].min(),'y_end':line_df['y'].max(),},ignore_index=True)

p.to_csv(f'output/line_parameter.txt',sep='\t',header=True,index=False)

