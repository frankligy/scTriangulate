import sys
import os
import math
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import rankdata
from scipy.sparse import issparse
import multiprocessing as mp
import logging
import scanpy as sc
import anndata as ad
import subprocess

from .shapley import approximate_shapley_value
from .metrics import *




def rank_pruning(sctri,discard=None,scale_sccaf=True,layer=None,assess_raw=False):

    if assess_raw:
        # construct data for shapley
        sctri.run_single_key_assessment(key='raw',scale_sccaf=scale_sccaf,layer=layer)
        subprocess.run(['rm','-r','{}'.format(os.path.join(sctri.dir,'scTriangulate_local_mode_enrichr/'))])
        score_info = copy.deepcopy(sctri.score['raw'])
        score_info.pop('cluster_to_doublet',None)  # don't consider doublet score

        data_shape0 = len(score_info['cluster_to_reassign'])
        data_shape1 = len(score_info)
        data = np.empty(shape=[data_shape0,data_shape1])
        df_index = list(score_info['cluster_to_reassign'].keys())
        df_columns = list(score_info.keys())
        for m,(metric,info) in enumerate(score_info.items()):
            for c,(cluster,score) in enumerate(info.items()):
                data[c,m] = score


        # compute shapley
        appro_shapley = approximate_shapley_value(data)
        df = pd.DataFrame(data=data,index=df_index,columns=df_columns)
        df['appro_shapley'] = appro_shapley
        df.sort_values(by='appro_shapley',ascending=False,inplace=True)


        # concat clusters that have only one cells back to the df for completeness
        all_clusters = sctri.adata.obs['raw'].unique()
        ones_clusters = list(set(all_clusters).difference(set(df.index)))
        data_ones = np.zeros(shape=[len(ones_clusters),data_shape1+1])  # add one for appro_shapley column
        df_ones = pd.DataFrame(data=data_ones,index=ones_clusters,columns=df.columns)
        df_total = pd.concat([df,df_ones])
        # min-max scaler the appro_shapley
        max_shapley = df_total['appro_shapley'].max()
        min_shapley = df_total['appro_shapley'].min()
        scaled = [(item - min_shapley) / (max_shapley - min_shapley) for item in df_total['appro_shapley'].values]
        df_total['scaled_shapley'] = scaled

    else:
        df_total = pd.DataFrame(index=sctri.adata.obs['raw'].unique())
        
    # add cluster size
    raw_size_dict = sctri.adata.obs['raw'].value_counts().to_dict()
    df_total['raw_cluster_size'] = df_total.index.map(raw_size_dict).values
    # add original cluster
    ori_size_dict = {}
    for key,clusters in sctri.size_dict.items():
        for cluster,size in clusters.items():
            label = key + '@' + cluster
            size = size
            ori_size_dict[label] = size
    df_total['ori_cluster_size'] = df_total.index.map(ori_size_dict).values
    # add winning fraction
    df_total['win_fraction'] = df_total['raw_cluster_size'].values / df_total['ori_cluster_size'].values

    if discard is None:
        modified_obs = copy.deepcopy(sctri.adata.obs)
        modified_obs['pruned'] = modified_obs['raw'].values
        return modified_obs,df_total
    else:
        assert discard < df_total.shape[0]
        invalid = df_total.index.values[-discard:].tolist()
        sctri.add_to_invalid(invalid)
        modified_obs, invalid_now = reassign_pruning(sctri,abs_thresh=0,remove1=False)
        confidence_dict = pd.Series(index=df_total.index,data=df_total['scaled_shapley'].values)
        modified_obs['confidence'] = modified_obs['pruned'].map(confidence_dict).values
        return modified_obs,df_total


def reassign_pruning(sctri,abs_thresh=10,remove1=True,reference=None): #if you want to build viewer,set it to True

    adata = sctri.adata
    obs = adata.obs
    invalid = copy.deepcopy(sctri.invalid)
    size_dict = sctri.size_dict
    marker_genes = sctri.uns['marker_genes']
    query = sctri.query
    if reference is None:
        reference = sctri.reference
    else:
        reference = reference

    # add too small clusters to invaild list as well
    obs['ori'] = np.arange(obs.shape[0])     
    vc = obs['raw'].value_counts()
    for key_cluster in vc.index:
        size = size_dict[key_cluster.split('@')[0]][key_cluster.split('@')[1]]
        if vc[key_cluster] < abs_thresh:
            invalid.append(key_cluster)
        # elif vc[key_cluster] < 0.05 * size:
        #     invalid.append(key_cluster)
    
    invalid = list(set(invalid))

    if len(invalid) == 0:     # no invalid clusters, so no need to do reassign anymore
        return obs,invalid

    # seperate valid and invalid, only operate on invalid
    valid_obs = obs.loc[~obs['raw'].isin(invalid),:]
    invalid_obs = obs.loc[obs['raw'].isin(invalid),:]
  
    # for invalid ones, find nearest centroid in valid one
    ## get pool
    num = 30
    pool = []
    for key in query:
        marker = marker_genes[key]
        for i in range(marker.shape[0]):
            used_marker_genes = marker.iloc[i]['purify']
            pick = used_marker_genes[:num]  # if the list doesn't have more than 30 markers, it is oK, python will automatically choose all
            pool.extend(pick)
    pool = list(set(pool))
    adata_now = adata[:,pool].copy()

    if issparse(adata_now.X):
        adata_now.X = adata_now.X.toarray()


    ## mean-centered and divide the std of the data
    tmp = adata_now.X
    from sklearn.preprocessing import scale
    tmp_scaled = scale(tmp,axis=0)
    adata_now.X = tmp_scaled

    ## reducing dimension 
    from sklearn.decomposition import PCA
    n_components = 30
    reducer = PCA(n_components=n_components)
    scoring = reducer.fit_transform(X=adata_now.X) 
    adata_reduced = ad.AnnData(X=scoring,obs=adata_now.obs,var=pd.DataFrame(index=['PC{}'.format(str(i+1)) for i in range(n_components)]))

    ## have train and test
    adata_train = adata_reduced[valid_obs.index.tolist(),:]
    adata_test = adata_reduced[invalid_obs.index.tolist(),:]

    ## get X,y for training
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    adata_train.obs['raw'] = adata_train.obs['raw'].astype('category')
    X = np.empty([len(adata_train.obs['raw'].cat.categories),scoring.shape[1]])
    y = []
    for i,cluster in enumerate(adata_train.obs['raw'].cat.categories):
        bool_index = adata_train.obs['raw']==cluster
        centroid = np.mean(adata_train.X[bool_index,:],axis=0)
        X[i,:] = centroid
        y.append(cluster)
    y = le.fit_transform(y)

    ## training
    from sklearn.neighbors import KNeighborsClassifier
    n_neighbors = 10
    if X.shape[0] < n_neighbors:
        n_neighbors = X.shape[0]
    model = KNeighborsClassifier(n_neighbors=n_neighbors,weights='distance')
    model.fit(X,y)

    ## predict invalid ones
    X_test = adata_test.X
    pred = model.predict(X_test)  # (n_samples,)
    result = le.inverse_transform(pred)

    # start to reassemble
    adata_test.obs['pruned'] = result
    adata_train.obs['pruned'] = adata_train.obs['raw']
    modified_obs = pd.concat([adata_train.obs,adata_test.obs])
    modified_obs.sort_values(by='ori',inplace=True)


    # for plotting purpose, any cluster within a reference = 1 will be reassigned to most abundant cluster
    if remove1:
        bucket = []
        for ref,subset in modified_obs.groupby(by=reference):
            vc2 = subset['pruned'].value_counts()
            most_abundant_cluster = vc2.loc[vc2==vc2.max()].index[0]  # if multiple, just pick the first one
            exclude_clusters = vc2.loc[vc2==1].index
            for i in range(subset.shape[0]):
                if subset.iloc[i]['pruned'] in exclude_clusters:
                    subset.loc[:,'pruned'].iloc[i] = most_abundant_cluster   # caution that Settingwithcopy issue
            bucket.append(subset)
        modified_obs = pd.concat(bucket)
        modified_obs.sort_values(by='ori',inplace=True)

    return modified_obs,invalid
    



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
    abs_thresh = 10 if obs.shape[0] < 50000 else 30
    for cluster in overlap_clusters:
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
        elif proportion_to_ref < 0.05 and count_cluster > abs_thresh:  # not nearly included, its fraction to reference is low, but absolute count is decent, keep it
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

    



