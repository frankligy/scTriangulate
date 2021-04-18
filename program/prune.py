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
    if fraction_r < 0.1 and fraction_c < 0.1:
        result = False
    else:
        result = True
    # next, assess if one is nearly included in reference
    if fraction_c > 0.9:
        nearly = True
    else:
        nearly = False
    return fraction_r,fraction_c,result,nearly



def reference_pruning(adata,reference):
    obs = adata.obs
    obs['ori'] = np.arange(obs.shape[0])     # keep original index order in one column
    pruned_chunks = [] # store pruned chunk, one chunk menas one reference cluster
    for chunk in obs.groupby(by=reference):
        subset = chunk[1]
        vc = subset['engraft'].value_counts()
        overlap_clusters = vc.index
        mapping = {}
        for cluster in overlap_clusters:
            r = {reference:chunk[0]}
            c = {cluster.split('_')[0]:cluster.split('_')[1]}
            fraction_r,fraction_c,result,nearly = inclusiveness(obs,r,c)
            if not result:  # no inclusive, go back to reference annotation
                mapping[cluster] = reference + '_' + chunk[0]
            else:
                proportion = vc.loc[cluster] / vc.sum()
                if proportion < 0.1 and not nearly:  # only cover < 5% reference cluster and it is not nearly included
                    mapping[cluster] = reference + '_' + chunk[0]
                else:
                    mapping[cluster] = cluster

        subset['reassign'] = subset['engraft'].map(mapping).values

        # change to most abundant type if engraft only have 1
        vc2 = subset['reassign'].value_counts()
        most_abundant_cluster = vc2.loc[vc2==vc2.max()].index[0]  # if multiple, just pick the first one
        exclude_clusters = vc2.loc[vc2==1].index
        for i in range(subset.shape[0]):
            if subset.iloc[i]['reassign'] in exclude_clusters:
                subset.loc[:,'reassign'].iloc[i] = most_abundant_cluster   # caution that Settingwithcopy issue
        pruned_chunks.append(subset)
    modified_obs = pd.concat(pruned_chunks)
    modified_obs.sort_values(by='ori',inplace=True)
    adata.obs = modified_obs

    



