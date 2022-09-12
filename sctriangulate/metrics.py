import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import gseapy as gp
import logging
import math
import os

from .logger import logger_sctriangulate
from .preprocessing import *


def check_filter_single_cluster(adata,key):
    vc = adata.obs[key].value_counts()
    exclude_clusters= vc.loc[vc==1].index
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth,:]
    return adata_valid

def doublet_compute(adata,key):
    cluster_to_doublet = {}
    for cluster in adata.obs[key].astype('category').cat.categories:
        mean_score = adata[adata.obs[key]==cluster,:].obs['doublet_scores'].values.mean()
        cluster_to_doublet[cluster] = mean_score
    return cluster_to_doublet


def compute_combo_score(rank_uns,cluster):
    rank_names = rank_uns['names'][cluster]
    rank_lfc = rank_uns['logfoldchanges'][cluster]
    rank_pval = rank_uns['pvals'][cluster]
    df = pd.DataFrame({'names':rank_names,'lfc':rank_lfc,'pval':rank_pval})
    # filter out down-regulated genes
    df = df.loc[df['lfc'] > 0, :]
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    # the rank of each gene by lfc, the larger, the better, make argsort result reverse
    temp = np.flip(np.argsort(df['lfc'].values))
    ranks_lfc = np.empty_like(temp)
    ranks_lfc[temp] = np.arange(len(df['pval'].values))
    # the rank of each gene by pval, the smaller, the better
    temp = np.argsort(df['pval'].values)
    ranks_pval = np.empty_like(temp)
    ranks_pval[temp] = np.arange(len(df['pval'].values))
    # combo rank score
    temp = (ranks_lfc + ranks_pval) / 2
    df['rank_lfc'] = ranks_lfc
    df['rank_pval'] = ranks_pval
    df['combo'] = temp
    df.sort_values(by='combo', inplace=True)
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    # filter out the genes if pval > 0.05
    df = df.loc[df['pval']<0.05,:]
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    return df

def run_enrichr(gene_list,key,name,folder,species,criterion):
    # run enrichr
    artifact = read_artifact_genes(species,criterion=1).reset_index()
    artifact_dict = artifact.groupby(by='class')['genes'].apply(lambda x:x.tolist()).to_dict()
    enr2 = gp.enrichr(gene_list=gene_list,
                    description=name,
                    gene_sets=artifact_dict,
                    background=20000,
                    outdir=os.path.join(folder,'scTriangulate_local_mode_enrichr'),
                    cutoff=0.1, # adj-p for plotting
                    verbose=True)
    enrichr_result = enr2.results
    enrichr_dict = {}
    for metric in artifact_dict.keys():
        if enrichr_result.shape[0] == 0:  # no enrichment for any of the above terms
            enrichr_dict[metric] = 0
        else:
            try:
                enrichr_score = -math.log10(enrichr_result.loc[enrichr_result['Term']==metric,:]['Adjusted P-value'].to_list()[0])
            except IndexError:
                enrichr_dict[metric] = 0
            else:
                enrichr_dict[metric] = enrichr_score
    return enrichr_dict

def run_gsea(gene_list,key,name,folder,species,criterion):
    artifact = read_artifact_genes(species,criterion=1).reset_index()
    artifact_dict = artifact.groupby(by='class')['genes'].apply(lambda x:x.tolist()).to_dict()
    artifact_dict_keys = list(artifact_dict.keys())
    df = pd.DataFrame({0: gene_list, 1: 1/(np.arange(len(gene_list))+1)}) # col 1 is for descending rank of gene
    gsea_dict = {}
    try:
        pre_res = gp.prerank(rnk=df, gene_sets=artifact_dict,
                            permutation_num=100,
                            outdir=os.path.join(folder,'scTriangulate_local_mode_gsea/{}/{}'.format(key,name)),
                            min_size=1,
                            max_size=10000,
                            seed=6,
                            verbose=True)   # run this will cause artifact dict decreasing !! Caveats!!!
    except:  # no hit return, all metrics are zero
        for metric in artifact_dict_keys:
            gsea_dict[metric] = (0,0)   # first is nes, second is #hit
    else:
        gsea_result = pre_res.res2d
        metric_get = set(gsea_result.index.tolist())
        for metric in artifact_dict_keys:
            if metric in metric_get:
                gsea_score = gsea_result.loc[gsea_result.index==metric,:]['nes'].to_list()[0]
                gsea_hits = gsea_result.loc[gsea_result.index==metric,:]['matched_size'].to_list()[0]
                gsea_dict[metric] = (gsea_score, gsea_hits)
            else:  # not enriched
                gsea_dict[metric] = (0,0)
    return gsea_dict


def read_artifact_genes(species,criterion):
    '''
    criterion1: all will be artifact
    criterion2: all will be artifact except cellcycle
    criterion3: all will be artifact except cellcycle, ribosome
    criterion4: all will be artifact except cellcycle, ribosome, mitochondrial
    criterion5: all will be artifact except cellcycle, ribosome, mitochondrial, antisense
    criterion6: all will be artifact except cellcycle, ribosome, mitochondrial, antisense, predict_gene
    '''
    artifact = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'artifact_genes.txt'),sep='\t',index_col=0)
    artifact = artifact.loc[artifact['species']==species,:]
    if criterion == 1:
        artifact = artifact
    elif criterion == 2:
        artifact = artifact.loc[~(artifact['class']=='cellcycle'),:]
    elif criterion == 3:
        artifact = artifact.loc[~((artifact['class']=='ribosome')|(artifact['class']=='cellcycle')),:]
    elif criterion == 4:
        artifact = artifact.loc[~((artifact['class']=='ribosome')|(artifact['class']=='cellcylce')|(artifact['class']=='mitochondrial')),:]
    elif criterion == 5:
        artifact = artifact.loc[~((artifact['class']=='ribosome')|(artifact['class']=='cellcylce')|(artifact['class']=='mitochondrial')|(artifact['class']=='antisense')),:]
    elif criterion == 6:
        artifact = artifact.loc[~((artifact['class']=='ribosome')|(artifact['class']=='cellcylce')|(artifact['class']=='mitochondrial')|(artifact['class']=='antisense')|(artifact['class']=='predict_gene')),:]
    return artifact

    

def purify_gene(genelist,species,criterion):
    result = []
    artifact = read_artifact_genes(species,criterion)
    artifact_genes = set(artifact.index.to_list())
    for gene in genelist:
        if gene not in artifact_genes:
            result.append(gene)
    return result

def marker_gene(adata, key, species, criterion, folder):
    # delete previous rank_gene_gruops if present
    if adata.uns.get('rank_genes_groups') != None:
        del adata.uns['rank_genes_groups']
    # perform t-test
    sc.tl.rank_genes_groups(adata, key, method='t-test',n_genes=adata.shape[1])
    all_genes = adata.var_names.values  # ndarray, all the genes
    all_clusters = adata.obs[key].cat.categories  # pd.Index, all the clusters
    cluster2gene = dict()  # {'cluster1':[gene1,gene2..]}
    rank_uns = adata.uns['rank_genes_groups']
    pre_computed_dfs = []
    for cluster in all_clusters:
        cluster2gene[cluster] = []
        df = compute_combo_score(rank_uns, cluster)
        pre_computed_dfs.append(df)
    for gene in all_genes:            
        index_store = []
        for i,cluster in enumerate(all_clusters):
            df = pre_computed_dfs[i]
            # get the rank of the gene in each cluster
            try:
                index = np.nonzero(df['names'].values == gene)[0][0]  # the rank of this gene in each cluster
            except IndexError:
                index = len(all_genes)
            index_store.append(index)
        if np.all(np.array(index_store) == len(all_genes)):
            continue
        assign = all_clusters[np.argmin(np.array(index_store))]  # get argmin, take the corresponding cluster
        cluster2gene[assign].append((gene,np.min(index_store)))
    # sort the cluster2gene
    for key_,value in cluster2gene.items():
        gene = [item[0] for item in value]
        rank = [item[1] for item in value]
        temp = sorted(zip(gene,rank),key=lambda x:x[1])
        cluster2gene[key_] = [item[0] for item in temp]
    result = pd.Series(cluster2gene).to_frame()
    result.columns = ['whole_marker_genes']

    '''
    now the result is a dataframe
                whole_marker_genes
    cluster1        gene_list
    cluster2        gene_list
    '''

    # now let's perform enrichr and GSEA, and get puried marker gene
    col_enrichr = []
    col_gsea = []
    col_purify = []   # genelist that have artifact genes removed
    logger_sctriangulate.warning('Running GSEAPY for marker genes for {}, requires Internet connection (not an error, just reminder)'.format(key))
    for cluster in result.index:
        enrichr_dict = run_enrichr(result.loc[cluster,:].to_list()[0],key=key,name=cluster,folder=folder,species=species,criterion=criterion)  # [0] because it is a [[gene_list]],we only need [gene_list]
        gsea_dict = run_gsea(result.loc[cluster,:].to_list()[0],key=key,name=cluster,folder=folder,species=species,criterion=criterion)
        purified = purify_gene(result.loc[cluster,:].to_list()[0],species,criterion) # the [0] is explained last line
        col_enrichr.append(enrichr_dict)
        col_gsea.append(gsea_dict)
        col_purify.append(purified)

    result['enrichr'] = col_enrichr
    result['gsea'] = col_gsea
    result['purify'] = col_purify

    return result


def reassign_score(adata,key,marker,regress_size=False):
    # get gene pool, slice the adata
    num = 30
    pool = []
    for i in range(marker.shape[0]):
        marker_genes = marker.iloc[i]['purify']
        pick = marker_genes[:num]  # if the list doesn't have more than 30 markers, it is oK, python will automatically choose all
        pool.extend(pick)
    pool = list(set(pool))
    adata_now = adata[:,pool].copy()

    # mean-centered and divide the std of the data
    tmp = adata_now.X
    from sklearn.preprocessing import scale
    tmp_scaled = scale(tmp,axis=0)
    adata_now.X = tmp_scaled
    
    # reducing dimension 
    from sklearn.decomposition import PCA
    reducer = PCA(n_components=30)
    scoring = reducer.fit_transform(X=adata_now.X) 

    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    scoring_y = le.fit_transform(adata_now.obs[key].astype('str'))
    order = le.classes_

    # compute the centroid of each cluster
    X = np.empty([len(adata_now.obs[key].cat.categories),scoring.shape[1]])
    y = []
    for i,cluster in enumerate(adata_now.obs[key].cat.categories):
        bool_index = adata_now.obs[key]==cluster
        centroid = np.mean(scoring[bool_index,:],axis=0)
        X[i,:] = centroid
        y.append(cluster)
    y = le.fit_transform(y)


    # train a KNN classifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import confusion_matrix
    # if number of centroid(training data) < N_neighbors, will raise error, we hard code it to be 10
    n_neighbors = 10
    if X.shape[0] < n_neighbors:
        n_neighbors = X.shape[0]
    model = KNeighborsClassifier(n_neighbors=n_neighbors,weights='distance')
    model.fit(X,y)
    pred = model.predict(scoring)  # (n_samples,)
    mat = confusion_matrix(scoring_y,pred)
    confusion_reassign = pd.DataFrame(data=mat,index=order,columns=order)
    accuracy = []
    for i in range(mat.shape[0]):
        accuracy.append(mat[i,i]/np.sum(mat[i,:]))
    cluster_to_accuracy = {}
    for i,cluster in enumerate(order):
        cluster_to_accuracy[cluster] = accuracy[i]

    # whether to regress out the clutser size effect or not
    if regress_size:
        key_metric_dict = cluster_to_accuracy
        key_size_dict = get_size_in_metrics(adata.obs,key)
        df_inspect = pd.concat([pd.Series(key_metric_dict),pd.Series(key_size_dict)],axis=1) # index is cluster, col1 is metric, col2 is size
        cluster_to_accuracy = regress_size(df_inspect,regressor='GLM',to_dict=True)

    del adata_now

    return cluster_to_accuracy, confusion_reassign

'''below is the part for regression score'''
def background_normalizer(df,n_neighbors=10,scale=True):
    # df is a two column dataframe where first column is metric and second column is size
    from copy import deepcopy
    df = deepcopy(df)
    df['order'] = np.arange(df.shape[0])
    col = []
    for i in range(df.shape[0]):
        this_metric = df[0][i]
        distance_to_this = (df[0] - this_metric).abs()
        df_tmp = deepcopy(df)
        df_tmp['distance'] = distance_to_this.values
        df_tmp.sort_values(by='distance',inplace=True)
        neighbors_metric = df_tmp.iloc[:,0][:n_neighbors].values
        mean_ = neighbors_metric.mean()
        std_ = neighbors_metric.std()
        if scale:
            if std_ == 0:
                col.append(0)
            else:
                col.append((this_metric-mean_)/std_)
        else:
            col.append(this_metric-mean_)
    df['normalized'] = col
    return df

def regress_size(df_inspect,regressor='background_zscore',n_neighbors=10,to_dict=False):
    
    # df_inspect, index is cluster name, col1 is metric, col2 is size
    if regressor == 'background_zscore':
        df_now = background_normalizer(df_inspect,n_neighbors,True)
        residual = df_now['normalized'].values
        df_inspect[0] = residual
        normalized_metric_series = df_inspect[0]
    elif regressor == 'background_mean':
        df_now = background_normalizer(df_inspect,n_neighbors,False)
        residual = df_now['normalized'].values
        df_inspect[0] = residual
        normalized_metric_series = df_inspect[0]
    elif regressor == 'GLM':
        endog = df_inspect[0] # metric
        exog = df_inspect[1]  # size
        import statsmodels.api as sm
        exog = sm.add_constant(exog,prepend=True)
        model = sm.GLM(endog,exog,family=sm.families.Gaussian())
        res = model.fit()
        residual = res.resid_response
        normalized_metric_series = residual
    elif regressor == 'Huber':
        endog = df_inspect[0] # metric
        exog = df_inspect[1]  # size
        from sklearn.linear_model import HuberRegressor
        model = HuberRegressor().fit(exog.values.reshape(-1,1),endog.values)
        prediction = model.predict(exog.values.reshape(-1,1))
        residual = endog.values - prediction
        # outliers = model.outliers_
        df_inspect[0] = residual
        normalized_metric_series = df_inspect[0]
    elif regressor == 'RANSAC':
        endog = df_inspect[0] # metric
        exog = df_inspect[1]  # size
        from sklearn.linear_model import RANSACRegressor
        model = RANSACRegressor().fit(exog.values.reshape(-1,1),endog.values)
        prediction = model.predict(exog.values.reshape(-1,1))
        residual = endog.values - prediction
        #outliers = np.logical_not(model.inlier_mask_)
        df_inspect[0] = residual
        normalized_metric_series = df_inspect[0]
    elif regressor == 'TheilSen':
        endog = df_inspect[0] # metric
        exog = df_inspect[1]  # size
        from sklearn.linear_model import TheilSenRegressor
        model = TheilSenRegressor().fit(exog.values.reshape(-1,1),endog.values)
        prediction = model.predict(exog.values.reshape(-1,1))
        residual = endog.values - prediction
        df_inspect[0] = residual
        normalized_metric_series = df_inspect[0]
    if to_dict:
        normalized_metric_dict = normalized_metric_series.to_dict()
        final = normalized_metric_dict
    else:
        final = normalized_metric_series
    return final




def tf_idf_bare_compute(df,cluster):
    '''
    now the df contains all the gene for and an additional column for cluster
    '''
    # compute its tf_idf
    tmp1 = df.loc[df['cluster'] == cluster, :].loc[:,df.columns!='cluster'].values  # (n_cells,n_genes)
    tf = np.count_nonzero(tmp1,axis=0) / tmp1.shape[0]  # (n_genes,)
    tf = tf + 1e-5
    tmp2 = df.loc[:,df.columns!='cluster'].values
    df_ = np.count_nonzero(tmp2,axis=0) / tmp2.shape[0]  # (n_genes,)
    df_ = df_ + 1e-5
    idf = -np.log10(df_)
    tf_idf_ori = tf * idf  # (n_genes,)
    return tf_idf_ori

def single_size_query(obs,c):
    # c would be {gs:ERP4}
    key = list(c.keys())[0]
    cluster = list(c.values())[0]
    size = obs.loc[obs[key]==cluster,:].shape[0]
    return size

def get_size_in_metrics(obs,key):
    key_size_dict = {}  # {ERP1:54,ERP2:100....}
    for cluster in obs[key].unique():
        size = single_size_query(obs,{key:cluster})
        key_size_dict[cluster] = size
    return key_size_dict

def tf_idf10_for_cluster(adata,key,species,criterion,regress_size=False,layer=None):
    if layer is None:
        df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)  
    else:
        df = pd.DataFrame(data=make_sure_mat_dense(adata.layers[layer]), index=adata.obs_names, columns=adata.var_names) 
    df['cluster'] = adata.obs[key].astype('str').values
    cluster_to_tfidf10 = {} # store tfidf10 score
    cluster_to_exclusive = {}   # store exclusivly expressed genes
    for item in adata.obs[key].cat.categories:
        a = tf_idf_bare_compute(df,item)
        a_names = adata.var_names
        test = pd.Series(data=a, index=a_names)
        test.sort_values(ascending=False, inplace=True)
        # remove artifact genes
        artifact = read_artifact_genes(species,criterion)
        artifact_genes = set(artifact.index.to_list())
        test_pure = test.loc[~test.index.isin(artifact_genes)]
        result10 = test_pure.iloc[9] 
        cluster_to_tfidf10[item] = result10
        cluster_to_exclusive[item] = test_pure.to_dict()
    exclusive_genes = pd.Series(cluster_to_exclusive,name='genes')

    # whether to regress out the clutser size effect or not
    if regress_size:
        key_metric_dict = cluster_to_tfidf10
        key_size_dict = get_size_in_metrics(adata.obs,key)
        df_inspect = pd.concat([pd.Series(key_metric_dict),pd.Series(key_size_dict)],axis=1) # index is cluster, col1 is metric, col2 is size   
        cluster_to_tfidf10 = regress_size(df_inspect,regressor='GLM',to_dict=True)

    del df
    
    return cluster_to_tfidf10, exclusive_genes


def tf_idf5_for_cluster(adata,key,species,criterion,regress_size=False,layer=None):
    if layer is None:
        df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)  
    else:
        df = pd.DataFrame(data=make_sure_mat_dense(adata.layers[layer]), index=adata.obs_names, columns=adata.var_names) 
    df['cluster'] = adata.obs[key].astype('str').values
    cluster_to_tfidf5 = {} # store tfidf1 score
    for item in adata.obs[key].cat.categories:
        a = tf_idf_bare_compute(df,item)
        a_names = adata.var_names
        test = pd.Series(data=a, index=a_names)
        test.sort_values(ascending=False, inplace=True)
        # remove artifact genes
        artifact = read_artifact_genes(species,criterion)
        artifact_genes = set(artifact.index.to_list())
        test_pure = test.loc[~test.index.isin(artifact_genes)]
        result5 = test_pure.iloc[4] 
        cluster_to_tfidf5[item] = result5

    # whether to regress out the clutser size effect or not
    if regress_size:
        key_metric_dict = cluster_to_tfidf5
        key_size_dict = get_size_in_metrics(adata.obs,key)
        df_inspect = pd.concat([pd.Series(key_metric_dict),pd.Series(key_size_dict)],axis=1) # index is cluster, col1 is metric, col2 is size
        cluster_to_tfidf5 = regress_size(df_inspect,regressor='GLM',to_dict=True)

    del df

    return cluster_to_tfidf5

def tf_idf1_for_cluster(adata,key,species,criterion,regress_size=False,layer=None):
    if layer is None:
        df = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)  
    else:
        df = pd.DataFrame(data=make_sure_mat_dense(adata.layers[layer]), index=adata.obs_names, columns=adata.var_names) 
    df['cluster'] = adata.obs[key].astype('str').values
    cluster_to_tfidf1 = {} # store tfidf1 score
    for item in adata.obs[key].cat.categories:
        a = tf_idf_bare_compute(df,item)
        a_names = adata.var_names
        test = pd.Series(data=a, index=a_names)
        test.sort_values(ascending=False, inplace=True)
        # remove artifact genes
        artifact = read_artifact_genes(species,criterion)
        artifact_genes = set(artifact.index.to_list())
        test_pure = test.loc[~test.index.isin(artifact_genes)]
        result1 = test_pure.iloc[0] 
        cluster_to_tfidf1[item] = result1

    # whether to regress out the clutser size effect or not
    if regress_size:
        key_metric_dict = cluster_to_tfidf1
        key_size_dict = get_size_in_metrics(adata.obs,key)
        df_inspect = pd.concat([pd.Series(key_metric_dict),pd.Series(key_size_dict)],axis=1) # index is cluster, col1 is metric, col2 is size
        cluster_to_tfidf1 = regress_size(df_inspect,regressor='GLM',to_dict=True)

    del df

    return cluster_to_tfidf1




def SCCAF_score(adata, key, species, criterion, scale_sccaf,regress_size=False):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix
    # define X and Y and remove artifact genes in the first place
    artifact = read_artifact_genes(species,criterion)
    artifact_genes = set(artifact.index.to_list())
    X = adata[:,~adata.var_names.isin(artifact_genes)].X.copy()  # from ArrayView to ndarray
    Y = adata.obs[key].values

    # mean-centered and divide the std of the data, if too many cells (>50000), no scale, liblinear solver is robust to unscaled data
    if scale_sccaf:
        tmp = X
        from sklearn.preprocessing import scale
        tmp_scaled = scale(tmp,axis=0)
        X = tmp_scaled
     
    # label encoding Y to numerical values
    le = LabelEncoder()
    Y = le.fit_transform(Y)
    # stratified split to traing and test, train and test, then get confusion matrix
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=0)
    for train_index, test_index in sss.split(X, Y):
        X_train = X[train_index]
        Y_train = Y[train_index]
        X_test = X[test_index]
        Y_test = Y[test_index]
    model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=100000)
    model.fit(X_train, Y_train)
    result = model.predict(X_test)
    m = confusion_matrix(Y_test, result)
    confusion_sccaf = pd.DataFrame(data=m,index=le.classes_,columns=le.classes_)
    # derive cluster reliability from confusion matrix for each cluster
    numeric2reliable = []  # [0.4,0.5...] length is the number of clusters involved in self-projection
    for i in range(m.shape[0]):
        numeric2reliable.append(m[i, i] / m[i, :].sum())
    cluster_to_SCCAF = {}
    for i in range(len(numeric2reliable)):
        cluster_to_SCCAF[le.classes_[i]] = numeric2reliable[i]

    # whether to regress out the clustser size effect or not
    if regress_size:
        key_metric_dict = cluster_to_SCCAF
        key_size_dict = get_size_in_metrics(adata.obs,key)
        df_inspect = pd.concat([pd.Series(key_metric_dict),pd.Series(key_size_dict)],axis=1) # index is cluster, col1 is metric, col2 is size
        cluster_to_SCCAF = regress_size(df_inspect,regressor='GLM',to_dict=True)

    del X
    del X_train
    del X_test

    return cluster_to_SCCAF, confusion_sccaf


















