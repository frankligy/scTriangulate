import sys
import os
import math
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
import multiprocessing as mp
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
import scanpy as sc
import anndata as ad



def reassign_pruning(adata):
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

# the following function will be used for referene pruning
def inclusiveness(obs,r,c):
    # r is the name of reference cluster, c is the name of cluster that overlap with r in the form of a dict
    # for example, r is {gs:ERP4}, c is {leiden1:16}
    r_key = list(r.keys())[0]
    r_cluster = list(r.values())[0]
    c_key = list(c.keys())[0]
    c_cluster = list(c.values())[0]
    obs[r_key] = obs[r_key].astype('str')
    obs[c_key] = obs[c_key].astype('str')
    # build set
    r_set = set(obs.loc[obs[r_key]==r_cluster,:].index.to_list())
    c_set = set(obs.loc[obs[c_key]==c_cluster,:].index.to_list())
    rc_i = r_set.intersection(c_set)
    fraction_r = len(rc_i)/len(r_set)
    fraction_c = len(rc_i)/len(c_set)
    return fraction_r,fraction_c


def run_reference_pruning(chunk,reference,size_dict,obs):
    subset = chunk[1]
    vc = subset['raw'].value_counts()
    overlap_clusters = vc.index
    mapping = {}
    for cluster in overlap_clusters:
        logging.info('process {}: {} contains {}'.format(os.getpid(),chunk[0],overlap_clusters))
        r = {reference:chunk[0]}
        c = {cluster.split('@')[0]:cluster.split('@')[1]}
        fraction_r,fraction_c = inclusiveness(obs,r,c)  # two cluster inclusive
        count_cluster = vc.loc[cluster]
        proportion_to_ref = vc.loc[cluster] / vc.sum()  # won cluster and reference inclusive
        proportion_to_self = vc.loc[cluster] / size_dict[cluster.split('@')[0]][cluster.split('@')[1]]
        if proportion_to_self >= 0.6:
            mapping[cluster] = cluster # nearly included, no matter how small its fraction is to the reference, keep it
        elif proportion_to_ref >= 0.05:
            mapping[cluster] = cluster  # not nearly included, but its fraction to reference is decent, keep it
        elif proportion_to_ref < 0.05 and count_cluster > 30:  # not nearly included, its fraction to reference is low, but absolute count is decent, keep it
            mapping[cluster] = cluster
        else:     # other wise, go back to reference cluster type
            mapping[cluster] = reference + '@' + chunk[0]

    subset['pruned'] = subset['raw'].map(mapping).values

    # change to most abundant type if pruned only have 1 cells, just for downstream DE analysis
    vc2 = subset['pruned'].value_counts()
    most_abundant_cluster = vc2.loc[vc2==vc2.max()].index[0]  # if multiple, just pick the first one
    exclude_clusters = vc2.loc[vc2==1].index
    for i in range(subset.shape[0]):
        if subset.iloc[i]['pruned'] in exclude_clusters:
            subset.loc[:,'pruned'].iloc[i] = most_abundant_cluster   # caution that Settingwithcopy issue
    return subset


def reference_pruning(obs,reference,size_dict):
    obs['ori'] = np.arange(obs.shape[0])     # keep original index order in one column
    pruned_chunks = [] # store pruned chunk, one chunk menas one reference cluster
    chunk_list = list(obs.groupby(by=reference))
    cores1 = len(chunk_list)
    cores2 = mp.cpu_count()
    cores = min(cores1,cores2)
    pool = mp.Pool(processes=cores)
    logging.info('spawn {} sub process for pruning'.format(cores))
    r = [pool.apply_async(run_reference_pruning,args=(chunk,reference,size_dict,obs)) for chunk in chunk_list]
    pool.close()
    pool.join()
    pruned_chunks = [collect.get() for collect in r]
    modified_obs = pd.concat(pruned_chunks)
    modified_obs.sort_values(by='ori',inplace=True)
    return modified_obs

    




