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
from scipy.sparse import issparse,csr_matrix
import multiprocessing as mp

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
    vc = adata.obs[key].astype('str').value_counts()
    exclude_clusters= vc.loc[vc==1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth,:]
    return adata_valid



# to parallelize, define singular function
def each_key_program(key):
    print(key,os.getpid())
    adata_to_compute = check_filter_single_cluster(adata,key)  # be a view
    print('finished filtering')
    result = marker_gene(adata_to_compute,key=key)
    print('finished marker gene computing')
    cluster_to_accuracy = reassign_score(adata_to_compute,key,result)
    print('finished reassign score')
    cluster_to_tfidf = tf_idf_for_cluster(adata_to_compute,key)
    print('finished tfidf')
    cluster_to_SCCAF = SCCAF_score(adata_to_compute,key)
    print('finished SCCAF')
    col_reassign = adata.obs[key].astype('str').map(cluster_to_accuracy).fillna(0).values
    col_tfidf = adata.obs[key].astype('str').map(cluster_to_tfidf).fillna(0).values
    col_SCCAF = adata.obs[key].astype('str').map(cluster_to_SCCAF).fillna(0).values

    # add a doublet score to each cluster
    cluster_to_doublet = doublet_compute(adata_to_compute,key)

    to_json = [cluster_to_accuracy,cluster_to_tfidf,cluster_to_SCCAF,cluster_to_doublet]
    to_viewer = list(cluster_to_accuracy.keys())

    # some diagnose plot
    draw_enrich_plots(key)
    draw_umap(adata,key)
    print('finished diagnostic')

    # all the intermediate results needed to be returned
    collect = {'key':key,
               'col_reassign':col_reassign,
               'col_tfidf':col_tfidf,
               'col_SCCAF':col_SCCAF,
               'to_json':to_json,
               'to_viewer':to_viewer}
    return collect


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



# give an adata, have raw attribute (only need X attribute), several obs column corresponding to different sets of annotations,
# users supplied umap if preferred

adata = sc.read('./input.h5ad')
query = ['ch_hsc','ch_baso','ch_ly6d','ch_prime','ch_thymo','un_cd127','un_hsc','un_kit','un_multilin','un_prime','un_thymo']
reference = 'ch_prime'

# precomputing size
size_dict,size_list = get_size(adata.obs,query)
c,s = size_sort(size_list)

if issparse(adata.X):
    adata.X = adata.X.toarray()


# add a doublet column
counts_matrix = adata.X
scrub = scr.Scrublet(counts_matrix)
doublet_scores,predicted_doublets = scrub.scrub_doublets(min_counts=1,min_cells=1)
adata.obs['doublet_scores'] = doublet_scores

sc.pl.umap(adata,color=['doublet_scores'],cmap='YlOrRd')
plt.savefig('./scTriangulate_diagnose/doublet.png',bbox_inches='tight')
plt.close()
print('finished doublet check')

del counts_matrix
del scrub


# compute metrics and map to original adata
data_to_json = {}
data_to_viewer = {}

results = []
for key in query:
    adata_to_compute = check_filter_single_cluster(adata,key)  # be a view
    print('finished filtering')
    result = marker_gene(adata_to_compute,key=key)
    print('finished marker gene computing')
    cluster_to_accuracy = reassign_score(adata_to_compute,key,result)
    print('finished reassign score')
    cluster_to_tfidf = tf_idf_for_cluster(adata_to_compute,key)
    print('finished tfidf')
    cluster_to_SCCAF = SCCAF_score(adata_to_compute,key)
    print('finished SCCAF')
    col_reassign = adata.obs[key].astype('str').map(cluster_to_accuracy).fillna(0).values
    col_tfidf = adata.obs[key].astype('str').map(cluster_to_tfidf).fillna(0).values
    col_SCCAF = adata.obs[key].astype('str').map(cluster_to_SCCAF).fillna(0).values

    # add a doublet score to each cluster
    cluster_to_doublet = doublet_compute(adata_to_compute,key)

    to_json = [cluster_to_accuracy,cluster_to_tfidf,cluster_to_SCCAF,cluster_to_doublet]
    to_viewer = list(cluster_to_accuracy.keys())

    # some diagnose plot
    draw_enrich_plots(key)
    draw_umap(adata,key)
    print('finished diagnostic')

    # all the intermediate results needed to be returned
    collect = {'key':key,
               'col_reassign':col_reassign,
               'col_tfidf':col_tfidf,
               'col_SCCAF':col_SCCAF,
               'to_json':to_json,
               'to_viewer':to_viewer}
    results.append(collect)

for collect in results:
    key = collect['key']
    adata.obs['reassign${}'.format(key)] = collect['col_reassign']
    adata.obs['tfidf${}'.format(key)] = collect['col_tfidf']
    adata.obs['SCCAF${}'.format(key)] = collect['col_SCCAF']
    data_to_json[key] = collect['to_json']
    data_to_viewer[key] = collect['to_viewer']

with open('./scTriangulate_present/score.json','w') as f:
    json.dump(data_to_json,f)
with open('./scTriangulate_present/key_cluster.p','wb') as f:
    pickle.dump(data_to_viewer,f)

adata.X = csr_matrix(adata.X)
adata.write('./scTriangulate_result/after_metrics_computing.h5ad')
adata.obs.to_csv('./scTriangulate_result/check_metrics.txt',sep='\t')
print('finished metrics computing and diagnose plot generaton')
print('----------------------------')


adata = sc.read('./scTriangulate_result/after_metrics_computing.h5ad')


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
    practical_colname = [name + '$' + key for name in score_colname]
    data[i,:,:] = adata.obs[practical_colname].values

# parallelize
def run_shapley(data):
    with open('./scTriangulate_present/log_shapley_step_{}.txt'.format(os.getpid()),'a') as log:
        log.write('This core needs to process {} cells\n'.format(data.shape[1]))
        log.flush()
        final = []
        intermediate = []
        for i in range(data.shape[1]):
            if i % 500 == 0 or i==data.shape[1]-1:
                print('Cell{}'.format(i),file=log,flush=True)
            layer = data[:,i,:]
            result = []
            for j in range(layer.shape[0]):
                result.append(shapley_value(j,layer))
            cluster_row = adata.obs.iloc[i].loc[query].values
            to_take = which_to_take(result,query,reference,cluster_row,size_dict)   # which annotation this cell should adopt
            final.append(to_take)    
            intermediate.append(result)
    return final,intermediate

final = []
intermediate = []
cores = mp.cpu_count()
sub_datas = np.array_split(data,cores,axis=1)  # [sub_data,sub_data,....]
pool = mp.Pool(processes=cores)
r = pool.map_async(run_shapley,sub_datas)
pool.close()
pool.join()
results = r.get()  # [(final,intermediate), (), ()...]
for collect in results:
    final.extend(collect[0])
    intermediate.extend(collect[1])
adata.obs['final_annotation'] = final
decisions = list(zip(*intermediate))
for i,d in enumerate(decisions):
    adata.obs['{}_shapley'.format(query[i])] = d
print('finished shapley computing')

adata.write('./scTriangulate_present/just_after_shapley.h5ad')

# adata = sc.read('./just_after_shapley.h5ad')

# assign
# parallelize
def run_assign(obs):       
    with open('./scTriangulate_present/log_assign_step_{}.txt'.format(os.getpid()),'a') as log:
        log.write('This core needs to process {} cells\n'.format(obs.shape[0]))
        log.flush()
        assign = []
        for i in range(obs.shape[0]):
            if i % 500 == 0 or i==obs.shape[0]-1:
                print('cell{}'.format(i),file=log,flush=True)
            name = obs.iloc[i,:].loc['final_annotation']
            cluster = obs.iloc[i,:].loc[name]
            concat = name + '$' + cluster
            assign.append(concat)   
    obs['engraft'] = assign
    return obs

obs = adata.obs
obs_index = np.arange(obs.shape[0])  # [0,1,2,.....]
cores = mp.cpu_count()
sub_indices = np.array_split(obs_index,cores)  # indices for each chunk [(0,1,2...),(56,57,58...),(),....]
sub_obs = [obs.iloc[sub_index,:] for sub_index in sub_indices]  # [sub_df,sub_df,...]
pool = mp.Pool(processes=cores)
r = pool.map_async(run_assign,sub_obs)
pool.close()
pool.join()
results = r.get()  # [sub_obs,sub_obs...]
obs = pd.concat(results)
adata.obs = obs
print('finished engraft')

adata.write('./scTriangulate_present/just_after_engraft.h5ad')



# prune
obs = reference_pruning(adata.obs,reference,size_dict)
adata.obs = obs
print('finished pruning')

adata.write('./scTriangulate_present/just_after_prune.h5ad')

# prefix with reference cluster
col1 = adata.obs['reassign']
col2 = adata.obs[reference]
col = []
for i in range(len(col1)):
    concat = reference + '$' + col2[i] + '|' + col1[i]
    col.append(concat)
adata.obs['reassign_prefix'] = col
print('finished prefix')

# print out
adata.obs.to_csv('./scTriangulate_present/shapley_annotation.txt',sep='\t')
adata.write('./scTriangulate_present/after_shapley_to_cellxgene.h5ad')


print('finished print out')

# inspection (DE, small umap, seperate adata h5ad)
plot_DE_umap_save(adata,reference)

print('finished inspection')

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

print('finished viewer building')








