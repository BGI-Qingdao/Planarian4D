
import os
import sys

import cv2
import numpy as np

from collections import Counter
import scipy.sparse
import pandas as pd

sys.path.append('/dellfsqd2/ST_OCEAN/USER/hankai/software/SpatialTranscript/site-packages')
from stomics import image


class Mask:
    def __init__(self, matrix):
        if isinstance(matrix, str):
            matrix = np.loadtxt(matrix, dtype=np.uint16)
        self.matrix = matrix
    def binning(self, bin_size):
        mtx= scipy.sparse.csc_matrix(tissue_mask)
        mtx = mtx.tocoo()
        tmp = []
        for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
            tmp.append([x, y, int(mask)])
        triplet = pd.DataFrame(tmp, columns=['x', 'y', 'mask'])
    
        triplet['xbin'] = (triplet.x / bin_size).astype(int) * bin_size
        triplet['ybin'] = (triplet.y / bin_size).astype(int) * bin_size
        triplet['bin'] = triplet.xbin.astype(str) + '_' + triplet.ybin.astype(str)

        index = [(-i, x) for i, x in enumerate(triplet['bin'].unique())]
        index = pd.DataFrame(index, columns=['N', 'bin'])

        triplet = triplet.merge(index, how='left', on='bin')
    
        matrix = np.zeros(shape=self.matrix.shape, dtype=int)
        matrix[triplet['x'], triplet['y']] = triplet['N']
        return matrix

    def to_binary(self):
        matrix = self.matrix
        mask = np.isin(matrix, [0], invert=True)
        matrix[mask] = 1
        return matrix

    def subtract(self, mask, label_area_cutoff=0.3):
        values = np.unique(self.matrix)
        if len(values) == 2:
            mask = cv2.bitwise_and(self.matrix, mask)
        else:
            self_mask = np.isin(self.matrix, [0], invert=True)
            mask = cv2.bitwise_and(self_mask, mask)

            orig_counter = Counter(self.matrix.flatten())

            mask = np.isin(mask, [0])
            filter_part = self.matrix[mask]
            filter_counter = Counter(filter_part.flatten())

            filter_labels = []
            for label, pixels in filter_counter.items():
                ratio = pixels / orig_counter[label]
                if ratio < label_area_cutoff:
                    continue
                filter_labels.append(label)

            filter_labels = list(set(filter_labels))
            mask = np.isin(self.matrix, filter_labels, invert=True)

        matrix = self.matrix[mask]
        return matrix

    def filter(self, on=None, cutoff=None):
        """label mask method
        * on: filter by minimum value of the input matrix
        """
        from scipy import ndimage
        labels = np.unique(self.matrix)
        medians = ndimage.median(
                self.matrix, 
                labels=labels, 
                index=np.arange(1, len(labels) + 1)
                )
        filtered_labels = [labels[i] for i, m in enumerate(medians) 
                if m > cutoff]

        mask = np.isin(self.matrix, filtered_labels)
        matrix = self.matrix[mask]
        return matrix

    def filter_by_diameter(self, min_size=1):
        """label mask method
        * min_size: max circo radius
        """
        dist_map = cv2.distanceTransform(self.mask, cv2.DIST_L2, cv2.DIST_MASK_PRECISE)
        _, radius, _, center = cv2.minMaxLoc(dist_map)

        filtered_cells = []
        for r, c in zip(radius, center):
            if r < min_size:
                filtered_cells.append(self.mask[c])
        mask = np.isin(self.mask, filtered_cells)
        self.mask[mask] = 0


def filter_by_size(cell_mask, size=10):
    filter_cells = []
    orig_counter = Counter(cell_mask.flatten())
    for cell, pixels in orig_counter.items():
        if pixels >= size:
            continue
        filter_cells.append(cell)
    mask = np.isin(cell_mask, filter_cells)
    cell_mask[mask] = 0
    return cell_mask

def shield_non_tissue(tissue_mask, cell_mask, min_ps=10, edge_cutoff=0.3):
    
    orig_counter = Counter(cell_mask.flatten())

    tissue_mask = np.isin(tissue_mask, [0])
    filter_part = cell_mask[tissue_mask]
    filter_counter = Counter(filter_part.flatten())

    filter_cells = []
    for cell, pixels in filter_counter.items():
        ratio = pixels / orig_counter[cell]
        if ratio < edge_cutoff:
            continue
        filter_cells.append(cell)

    filter_cells = list(set(filter_cells))
    mask = np.isin(cell_mask, filter_cells)
    cell_mask[mask] = 0
    return cell_mask

def filter_blur_cell(blur_mask, cell_mask, cutoff=30):
    blur_mask = blur_mask < 30
    cv2.imwrite('blur_mask.png', blur_mask.astype(int))
    #np.savetxt('blur.txt', blur_mask, fmt='%d')
    cell_mask[blur_mask] = 0
    return cell_mask

def complete_non_cell(tissue_mask, label_cell_mask, bin_size=20, edge_cutoff=0.3):
    
    mtx= scipy.sparse.csc_matrix(tissue_mask)
    mtx = mtx.tocoo()
    tmp = []
    for x, y, mask in zip(mtx.row, mtx.col, mtx.data):
        tmp.append([x, y, int(mask)])
    triplet = pd.DataFrame(tmp, columns=['x', 'y', 'mask'])
    
    triplet['xbin'] = (triplet.x / bin_size).astype(int) * bin_size
    triplet['ybin'] = (triplet.y / bin_size).astype(int) * bin_size
    triplet['bin'] = triplet.xbin.astype(str) + '_' + triplet.ybin.astype(str)

    index = [(-i, x) for i, x in enumerate(triplet['bin'].unique())]
    index = pd.DataFrame(index, columns=['N', 'bin'])

    triplet = triplet.merge(index, how='left', on='bin')
    
    mix_mask = np.zeros(shape=tissue_mask.shape, dtype=int)
    mix_mask[triplet['x'], triplet['y']] = triplet['N']
    mix_counter = Counter(mix_mask.flatten())
    
    cell_mask = np.isin(label_cell_mask, [0], invert=True)

    overlap_mix = mix_mask[cell_mask]
    overlap_counter = Counter(overlap_mix.flatten())

    filter_bin = []
    for bin_name, pixels in overlap_counter.items():
        ratio = pixels / mix_counter[bin_name]
        if ratio < edge_cutoff:
            continue
        filter_bin.append(bin_name)
    mask = np.isin(mix_mask, filter_bin)
    mix_mask[mask] = 0

    #mix_mask[cell_mask] = label_cell_mask
    np.copyto(mix_mask, label_cell_mask, casting='unsafe', where=cell_mask)

    return mix_mask

if __name__ == '__main__':

    tissue_mask = cv2.imread(sys.argv[1], cv2.IMREAD_UNCHANGED)
    cell_mask = np.loadtxt(sys.argv[2])
    blur_mask = np.loadtxt(sys.argv[3])
    img = cv2.imread(sys.argv[4])

    cell_mask = filter_by_size(cell_mask, size=10)
    _, fimg = image.overlayoutlines(img, cell_mask, draw=True)
    cv2.imwrite('filter_cell.png', fimg)

    cell_mask = shield_non_tissue(tissue_mask, cell_mask)
    _, simg = image.overlayoutlines(img, cell_mask, draw=True)
    cv2.imwrite('shield_cell.png', simg)

    cell_mask = filter_blur_cell(blur_mask, cell_mask)
    _, bimg = image.overlayoutlines(img, cell_mask, draw=True)
    cv2.imwrite('blur_cell.png', bimg)

    mix_mask = complete_non_cell(tissue_mask, cell_mask)
    np.savetxt('mix_mask.txt', mix_mask, fmt='%d')
    _, mimg = image.overlayoutlines(img, mix_mask, draw=True)
    cv2.imwrite('mix_mask.png', mimg)


