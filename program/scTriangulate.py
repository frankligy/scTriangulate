#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import sys
import os
import math
import json
import pickle
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata

import scanpy as sc
import anndata as ad
import gseapy as gp
import scrublet as scr

from shapley import *
from metrics import *
from diagnose import *
from viewer import *
from prune import *
from inspection import *


# helper functions for running main program
def check_filter_single_cluster(adata,key):
    vc = adata.obs['gs'].astype('str').value_counts()
    exclude_clusters= vc.loc[vc==1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth,:]
    return adata_valid



# set some file path
if not os.path.exists('./scTriangulate_result'):
    os.mkdir('./scTriangulate_result')
if not os.path.exists('./scTriangulate_diagnose'):
    os.mkdir('./scTriangulate_diagnose')
if not os.path.exists('./scTriangulate_local_mode_enrichr'):
    os.mkdir('./scTriangulate_local_mode_enrichr')
if not os.path.exists('./scTriangulate_present'):
    os.mkdir('./scTriangulate_present')
if not os.path.exists('./scTriangulate_inspection'):
    os.mkdir('./scTriangulate_inspection')



# give an adata, have raw attribute, several obs column corresponding to different sets of annotations

adata = sc.read('./leiden_gs_nathan.h5ad')
query = ['leiden0.5','leiden1','leiden2','gs']
reference = 'gs'

# add a doublet column
counts_matrix = adata.X.copy()
scrub = scr.Scrublet(counts_matrix)
doublet_scores,predicted_doublets = scrub.scrub_doublets(min_counts=1,min_cells=1)
adata.obs['doublet_scores'] = doublet_scores
print('finished doublet check')

# compute metrics and map to original adata
data_to_json = {}
data_to_viewer = {}
for key in query:
    print(key)
    adata_to_compute = check_filter_single_cluster(adata,key)  # every adata will be a copy
    print('finished filtering')
    result = marker_gene(adata_to_compute,key=key)
    print('finished marker gene computing')
    result.to_csv('./scTriangulate_result/marker_{0}.txt'.format(key),sep='\t')
    cluster_to_accuracy = reassign_score(adata_to_compute,key,result)
    print('finished reassign score')
    cluster_to_tfidf = tf_idf_for_cluster(adata_to_compute,key)
    print('finished tfidf')
    cluster_to_SCCAF = SCCAF_score(adata_to_compute,key)
    print('finished SCCAF')
    adata.obs['reassign_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_accuracy).fillna(0).values
    adata.obs['tfidf_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_tfidf).fillna(0).values
    adata.obs['SCCAF_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_SCCAF).fillna(0).values

    # add a doublet score to each cluster
    cluster_to_doublet = doublet_compute(adata_to_compute,key)

    data_to_json[key] = [cluster_to_accuracy,cluster_to_tfidf,cluster_to_SCCAF,cluster_to_doublet]
    data_to_viewer[key] = list(cluster_to_accuracy.keys())

with open('./scTriangulate_present/score.json','w') as f:
    json.dump(data_to_json,f)

with open('./scTriangulate_present/key_cluster.p','wb') as f:
    pickle.dump(data_to_viewer,f)

adata.write('./scTriangulate_result/after_metrics_computing.h5ad')
adata.obs.to_csv('./scTriangulate_result/check_metrics.txt',sep='\t')
print('finished metrics computing')
print('----------------------------')

#adata = sc.read('./scTriangulate_result/after_metrics_computing.h5ad')

# draw diagnostic plots
sc.pl.umap(adata,color=['doublet_scores'],cmap='YlOrRd')
plt.savefig('./scTriangulate_diagnose/doublet.png',bbox_inches='tight')
plt.close()

for key in query:
    draw_enrich_plots(key)
    draw_umap(adata,key)

print('finished diagnose')


# compute shaley value
score_colname = ['reassign','tfidf','SCCAF']
data = np.empty([len(query),adata.obs.shape[0],len(score_colname)])  # store the metric data for each cell
'''
data:
depth is how many sets of annotations
height is how many cells
width is how many score metrics
'''
for i,key in enumerate(query):
    practical_colname = [name + '_' + key for name in score_colname]
    data[i,:,:] = adata.obs[practical_colname].values

final = []
intermediate = []
for i in range(data.shape[1]):
    layer = data[:,i,:]
    result = []
    for j in range(layer.shape[0]):
        result.append(shapley_value(j,layer))
    to_take = which_to_take(adata,result,query,reference)   # which annotation this cell should adopt
    final.append(to_take)    
    intermediate.append(result)
adata.obs['final_annotation'] = final
decisions = zip(*intermediate)
for i,d in enumerate(decisions):
    adata.obs['{}_shapley'.format(query[i])] = d

print('finished shapley computing')

# assign
assign = []
for i in range(adata.obs.shape[0]):
    name = adata.obs.iloc[i,:].loc['final_annotation']
    cluster = adata.obs.iloc[i,:].loc[name]
    concat = name + '_' + cluster
    assign.append(concat)   
adata.obs['engraft'] = assign

print('finished engraft')

# prune
reference_pruning(adata,reference)
print('finished pruning')

# prefix with reference cluster
col1 = adata.obs['reassign']
col2 = adata.obs[reference]
col = []
for i in range(len(col1)):
    concat = reference + '_' + col2[i] + '|' + col1[i]
    col.append(concat)
adata.obs['reassign_prefix'] = col

# print out
adata.obs.to_csv('./scTriangulate_present/shapley_annotation.txt',sep='\t')
adata.write('./scTriangulate_present/after_shapley.h5ad')
adata.raw.to_adata().write('./scTriangulate_present/after_shapley_to_cellxgene.h5ad')

print('finished print out')

# inspection (DE, small umap, seperate adata h5ad)
adata = sc.read('./scTriangulate_present/after_shapley.h5ad')
plot_DE_umap_save(adata,reference)

# plot
fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
sc.pl.umap(adata,color='final_annotation',frameon=False,ax=ax[0])
sc.pl.umap(adata,color='final_annotation',frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
plt.savefig('./scTriangulate_present/umap_shapley_final_annotation.pdf',bbox_inches='tight')
plt.close()


fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})   # for reassign
sc.pl.umap(adata,color='reassign',frameon=False,ax=ax[0])
sc.pl.umap(adata,color='reassign',frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
plt.savefig('./scTriangulate_present/umap_shapley_reassign.pdf',bbox_inches='tight')
plt.close()

print('finished plotting')

# scTriangulate viewer
with open('./scTriangulate_present/key_cluster.p','rb') as f:
    data_to_viewer = pickle.load(f)
with open('./scTriangulate_present/score.json','r') as f:
    data_to_json = json.load(f)
with open('./scTriangulate_diagnose/viewer.html','w') as f:
    f.write(to_html(data_to_viewer,data_to_json))
with open('./scTriangulate_inspection/inspection.html','w') as f:
    f.write(inspection_html(data_to_viewer,reference))

os.system('cp ./viewer/viewer.css ./scTriangulate_diagnose')
os.system('cp ./viewer/viewer.js ./scTriangulate_diagnose')
os.system('cp ./viewer/inspection.css ./scTriangulate_inspection')
os.system('cp ./viewer/inspection.js ./scTriangulate_inspection')








