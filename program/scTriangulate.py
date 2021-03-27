#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import gseapy as gp
import math
from metrics import *
from diagnose import *

def check_filter_single_cluster(adata,key):
    vc = adata.obs['gs'].astype('str').value_counts()
    exclude_clusters= vc.loc[vc==1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth,:]
    return adata_valid

def shapley_value(index,data):
    '''
    it takes in a 2d ndarray, row means how many players, col means the metrics for a game performance
    [[0.5,0.4,0.7],
     [2.3,4.5,6.7],
     [0.1,0.1,0.1],
     [9.1,9.2,9.9]]

    also need to take in an index, row index, tell function which player's contribution you want to compute
    '''
    from itertools import combinations
    import math
    others = np.delete(data,index,axis=0)  # others except the one that we want to query
    others_index_to_combine = np.arange(others.shape[0])  # others' index in others matrix
    all_combine = []   # [(0,),(1,),(0,1)...]  # all the combinations of others' coalitions
    for r in range(others.shape[0]):
        coalition = list(combinations(others_index_to_combine,r=r+1))
        all_combine.extend(coalition)
    shapley = 0
    for item in all_combine:  # loop on each coalitions, compute normalized surplus
        feature_dim = data.shape[1]

        others_matrix = others[np.array(item),:].reshape(-1,feature_dim)
        player_row = data[index,:].reshape(-1,feature_dim)
        total_matrix = np.concatenate([others_matrix,player_row],axis=0)
        index_matrix = np.argsort(total_matrix,axis=0)
        player_rank = index_matrix[-1,:]
        good_or_not = (player_rank == total_matrix.shape[0]-1)
        player_rank[np.logical_not(good_or_not)] = 0
        surplus = player_rank.sum()

        s = len(item)
        n = data.shape[0]
        normalize_constant = math.factorial(s) * math.factorial(n-s-1) / math.factorial(n)
        value = normalize_constant * surplus
        shapley += value
    return shapley


def pruning(adata):
    # recompute the distances as a safe bet
    sc.pp.neighbors(adata,n_neighbors=10,n_pcs=30)
    distances = adata.obsp['distances'].toarray()
    # label red flag cluster
    vc = adata.obs['engraft'].astype('str').value_counts()
    red_flag = vc[vc < 6].index.values   # a ndarray storing cluster names that will be removed
    # traverse cells
    pending_modify = adata.obs['engraft'].astype('str').values
    for i in range(adata.obs.shape[0]):
        cluster = adata.obs['engraft'].astype('str').iloc[i]
        if cluster in red_flag:  # need to reassign
            neighbors_indices = np.nonzero(distances[i,:])[0]   # K nearest neighbors' indices
            neighbors_clusters = adata.obs['engraft'].astype('str').iloc[neighbors_indices].values  # the cluster name of all the neighbors
            # remove neighbors that actually are red flag as well
            neighbors_clusters = neighbors_clusters[[False if nbr in red_flag else True for nbr in neighbors_clusters]]
            u,counts = np.unique(neighbors_clusters,return_counts=True)
            assignment = sorted(zip(u,counts),key=lambda x:x[1],reverse=True)[0][0]  # assign to the most frequence cluster among its neighbors
            pending_modify[i] = assignment
            #print(cluster,neighbors_clusters,assignment)
        else:
            continue
    # add back
    adata.obs['reassign'] = pending_modify


# give an adata, have raw attribute, several obs column corresponding to different sets of annotations

adata = sc.read('./leiden_gs.h5ad')
query = ['leiden0.5','leiden1','leiden2','gs']



# compute metrics and map to original adata
data_to_json = {}
for key in query:
    print(key)
    adata_to_compute = check_filter_single_cluster(adata,key)  # every adata will be a copy
    print('finished filtering')
    result = marker_gene(adata_to_compute,key=key)
    print('finished marker gene computing')
    result.to_csv('./marker_{0}.txt'.format(key),sep='\t')
    cluster_to_accuracy = reassign_score(adata_to_compute,key,result)
    print('finished reassign score')
    cluster_to_tfidf = tf_idf_for_cluster(adata_to_compute,key)
    print('finished tfidf')
    cluster_to_SCCAF = SCCAF_score(adata_to_compute,key)
    print('finished SCCAF')
    adata.obs['reassign_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_accuracy).fillna(0).values
    adata.obs['tfidf_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_tfidf).fillna(0).values
    adata.obs['SCCAF_{}'.format(key)] = adata.obs[key].astype('str').map(cluster_to_SCCAF).fillna(0).values

    data_to_json[key] = [cluster_to_accuracy,cluster_to_tfidf,cluster_to_SCCAF]

with open('./score.json','w') as f:
    json.dump(data_to_json,f)

adata.write('./after_metrics_computing.h5ad')
adata.obs.to_csv('check_metrics.txt',sep='\t')
print('finished metrics computing')


# draw diagnostic plots
'''
using marker gene and exclusive gene result ouput
generate a folder where an index.html and cluster_1.html, cluster_proNeu_1.html
each html will have a barplot showing artifact gene enrichment, and UMAP for top10 marker gene and top10 exclusive genes
'''
for key in query:
    draw_enrich_plots(key)
    draw_umap(adata,key)

print('finished diagnose')
sys.exit('stop')

# compute shaley value
adata = sc.read('./after_metrics_computing.h5ad')
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
    to_take = query[result.index(max(result))]   # which annotation this cell should adopt
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

print('finished reassign')

# prune
pruning(adata)

print('finished pruning')

# print out
adata.obs.to_csv('./shapley_annotation.txt',sep='\t')
adata.write('./after_shapley.h5ad')

print('finished print out')

# plot
fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
sc.pl.umap(adata,color='final_annotation',frameon=False,ax=ax[0])
sc.pl.umap(adata,color='final_annotation',frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
plt.savefig('./umap_shapley_final_annotation.pdf',bbox_inches='tight')
plt.close()


fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})   # for reassign
sc.pl.umap(adata,color='reassign',frameon=False,ax=ax[0])
sc.pl.umap(adata,color='reassign',frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
plt.savefig('./umap_shapley_reassign.pdf',bbox_inches='tight')
plt.close()

print('finished plotting')






