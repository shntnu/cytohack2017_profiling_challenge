#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:51:34 2017

@author: ajaver
"""
#a) When our methods tell us a profile differs from the rest, how can we better 
#explore the causes of this difference, aside from just starting at rank=ordered
#list of features

#b) How can we better classify conditions in an unsupervised manner

#c) can we related chemical data to morphological data

import pandas as pd
import os

main_dir = '/Users/ajaver/Downloads/bbbc022'

fnames = [x for x in os.listdir(main_dir) if x.endswith('.csv')]


all_data  = []
for fname in fnames:
    df = pd.read_csv(os.path.join(main_dir, fname))
    all_data.append(df)
#%%
import numpy as np
all_data_d = pd.concat(all_data)
all_data_d.index = np.arange(all_data_d.shape[0])
#%%

feat_columns = [x for x in all_data_d.columns if 'Metadata' not in x]
#%%

all_data_n = all_data_d.copy()
for plate_id, data_plate in all_data_d.groupby('Metadata_Plate'):
    print(plate_id)
    control_columns = data_plate['Metadata_pert_type'] == 'control'
    controls_feats = data_plate.loc[control_columns.values, feat_columns]
    control_mean = controls_feats.mean()
    control_std = controls_feats.std()
    all_data_n.loc[data_plate.index, feat_columns] = (data_plate[feat_columns]-control_mean)/control_std
#%%
#remove features with nans (the lowest has 2811)
n_bads = all_data_n[feat_columns].isnull().sum()
bad_feats= n_bads[n_bads>0].index
good_cols = [x for x in all_data_n.columns if not x in bad_feats]
good_feats = [x for x in good_cols if x in feat_columns]

all_data_n = all_data_n[good_cols]
all_data_n.to_csv('bbbc022_normalized.csv')
#%%
good_data = all_data_n['Metadata_mmoles_per_liter'].isin((0,5)).sum()
all_data_n['Metadata_compound_name'] = all_data_n['Metadata_compound_name'].str.lower()
for compound, compound_data in all_data_n.groupby('Metadata_compound_name'):
    compound
#%%
#all_data_n = all_data_n[bad_feats]
#%%
import matplotlib.pylab as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

X = all_data_n[good_feats]

pca = PCA(n_components=20)
pca.fit(X)
Y_pca = pca.transform(all_data_n[good_feats])
plt.plot(np.cumsum(pca.explained_variance_ratio_))
#%%
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
#%%
plt.figure()

tsne = TSNE(n_components=2)
Y = tsne.fit_transform(Y_pca) 
#%%
#tsne = TSNE(n_components=2)
#Y2 = tsne.fit_transform(X) 
#%%
#all_data_n.to_csv('bbbc022_normalized.csv')
#%%
all_data_p = pd.read_csv('bbbc022_normalized_prints.csv')

#%%
all_data_r = all_data_p[~all_data_p['Metadata_smiles_fingerprint'].isnull()]
finger_d = all_data_r['Metadata_smiles_fingerprint']
#%%
#%%
def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

dd = min([hamming_distance(finger_d[0], X) for X in finger_d])
#%%
import seaborn

#%%
all_rows = []
for ii, row in all_data_p.iterrows():
    mm = row['Metadata_compound_name']
    row_name = mm.lower() if isinstance(mm, str) else 'control' + '_{:.3}'.format(row['Metadata_mmoles_per_liter'])
    for feat in good_feats:
        dd = (row_name, feat, row[feat])
        all_rows.append(dd)
    


#%%

new_map = []
#%%
for x in all_data_p.columns:
    print(x)

#%%
#get_centroid = lambda x : [str(int(x)) for x in np.round(np.random.rand(2048))]
#centroids = [get_centroid(x) for x in range(10)]
#
#data = finger_d
#k = len(centroids)
#while True:
#    distances = [[hamming_distance(x,y) for x in centroids]for y in data]
#    data_to_centroids = [min(enumerate(x),key=lambda x:x[1])[0] for x in distances]
#    break
##%%
#    dd =  [
#				map(
#					lambda x: data[x[0]], # get the actual data point
#					filter(
#						lambda x: x[1]==y,
#						enumerate(data_to_centroids)
#					) # get the points of cluster y as tuples (data_idx,cluster)
#				)
#				for y in range(k) # for each temporary cluster
#			]


#from kmeans import kmeans_hamming
#isinstance(x, str) else x for x in all_data_p['Metadata_smiles_fingerprint']




#%%
'''
Metadata_Plate                                                   20585
Metadata_Well                                                      A01
Metadata_Assay_Plate_Barcode                                     20585
Metadata_Plate_Map_Name                                   H-BIOA-002-1
Metadata_well_position                                             A01
Metadata_broad_sample                           BRD-K98763141-001-06-8
Metadata_mmoles_per_liter                                       8.8584
Metadata_source_name                         Biomol International Inc.
Metadata_compound_name                                   NIFLUMIC ACID
Metadata_smiles                      OC(=O)c1cccnc1Nc2cccc(c2)C(F)(F)F
Metadata_solvent                                                  DMSO
Metadata_pert_id                                         BRD-K98763141
Metadata_pert_mfc_id                            BRD-K98763141-001-06-8
Metadata_pert_well                                                 A01
Metadata_pert_id_vendor                                            NaN
Metadata_cell_id                                               unknown
Metadata_broad_sample_type                                         trt
Metadata_pert_vehicle                                             DMSO
Metadata_pert_type                                                 trt

'''