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

def run_gsea_and_enrichr(gene_list,name):
    # run enrichr
    artifact = pd.read_csv('./artifact_genes.txt',sep='\t')
    artifact_dict = artifact.groupby(by='notes')['genes'].apply(lambda x:x.tolist()).to_dict()
    enr2 = gp.enrichr(gene_list=gene_list,
                    description=name,
                    gene_sets=artifact_dict,
                    background=20000,
                    outdir='./scTriangulate_local_mode_enrichr/',
                    cutoff=0.1, # adj-p for plotting
                    verbose=True)
    enrichr_result = enr2.results
    enrichr_dict = {}
    for metric in ['Sex_associate', 'mitochondrial', 'predict_gene_Gm', 'predict_gene_dot', 'ribosomal_Rpl', 'ribosomal_Rps']:
        if enrichr_result.shape[0] == 0:  # no enrichment for any of the above terms
            enrichr_dict[metric] = 0
        else:
            try:
                enrichr_score = -math.log10(enrichr_result.loc[enrichr_result['Term']==metric,:]['Adjusted P-value'].to_list()[0])
            except IndexError:
                enrichr_dict[metric] = 0
            else:
                enrichr_dict[metric] = enrichr_score
        
    # run GSEA
    artifact = pd.read_csv('./artifact_genes.txt',sep='\t')
    artifact_dict = artifact.groupby(by='notes')['genes'].apply(lambda x:x.tolist()).to_dict()
    df = pd.DataFrame({0: gene_list, 1: 1/(np.arange(len(gene_list))+1) })
    gsea_dict = {}
    try:
        pre_res = gp.prerank(rnk=df, gene_sets=artifact_dict,
                            permutation_num=100,
                            outdir='./scTriangulate_local_mode_prerank/{}'.format(name),
                            min_size=1,
                            max_size=10000,
                            seed=6,
                            verbose=True)
    except:  # no hit return, all metrics are zero
        for metric in ['Sex_associate', 'mitochondira', 'predict_gene_Gm', 'predict_gene_dot', 'ribosomal_Rpl', 'ribosomal_Rps']:
            gsea_dict[metric] = 0
    else:
        gsea_result = pd.read_csv('./scTriangulate_local_mode_prerank/{}/gseapy.prerank.gene_sets.report.csv'.format(name))
        for metric in ['Sex_associate', 'mitochondira', 'predict_gene_Gm', 'predict_gene_dot', 'ribosomal_Rpl', 'ribosomal_Rps']:
            try:
                gsea_score = -math.log10(gsea_result.loc[gsea_result['Term']==metric,:]['pval'].to_list()[0])
            except:   # some pval=0, and pval may not be a good way to reflect that.
                gsea_dict[metric] = 0
            else:
                gsea_dict[metric] = gsea_score
    return enrichr_dict,gsea_dict

def purify_gene(genelist):
    result = []
    artifact = pd.read_csv('./artifact_genes.txt',sep='\t')
    artifact_genes = set(artifact['genes'].to_list())
    for gene in genelist:
        if gene not in artifact_genes:
            result.append(gene)
    return result

def marker_gene(adata, key):
    # delete previous rank_gene_gruops if present
    if adata.uns.get('rank_genes_groups') != None:
        del adata.uns['rank_genes_groups']
    # perform t-test
    sc.tl.rank_genes_groups(adata, key, method='t-test',n_genes=adata.raw.shape[1])
    all_genes = adata.raw.var_names.values  # ndarray, all the genes
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
    for key,value in cluster2gene.items():
        gene = [item[0] for item in value]
        rank = [item[1] for item in value]
        temp = sorted(zip(gene,rank),key=lambda x:x[1])
        cluster2gene[key] = [item[0] for item in temp]
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
    for cluster in result.index:
        enrichr_dict, gsea_dict = run_gsea_and_enrichr(result.loc[cluster,:].to_list()[0],name=cluster)  # [0] because it is a [[gene_list]],we only need [gene_list]
        purified = purify_gene(result.loc[cluster,:].to_list()[0]) # the [0] is explained last line
        col_enrichr.append(enrichr_dict)
        col_gsea.append(gsea_dict)
        col_purify.append(purified)

    result['enrichr'] = col_enrichr
    result['gsea'] = col_gsea
    result['purify'] = col_purify

    return result


def reassign_score(adata,key,marker):
    # get gene pool, slice the adata
    num = 30
    pool = []
    for i in range(marker.shape[0]):
        marker_genes = marker.iloc[i]['purify']
        pick = marker_genes[:num]  # if the list doesn't have more than 30 markers, it is oK, python will automatically choose all
        pool.extend(pick)
    pool = list(set(pool))
    adata_now = adata.raw.to_adata()
    adata_now = adata_now[:,pool]
    
    # reducing dimension 
    from sklearn.decomposition import PCA
    reducer = PCA(n_components=30)
    scoring = reducer.fit_transform(X=adata_now.X.toarray())

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
    model = KNeighborsClassifier(n_neighbors=10,weights='distance')
    model.fit(X,y)
    pred = model.predict(scoring)  # (n_samples,)
    mat = confusion_matrix(scoring_y,pred)
    pd.DataFrame(data=mat,index=order,columns=order).to_csv('./scTriangulate_result/confusion_reassign_{}.txt'.format(key),sep='\t')
    accuracy = []
    for i in range(mat.shape[0]):
        accuracy.append(mat[i,i]/np.sum(mat[i,:]))
    cluster_to_accuracy = {}
    for i,cluster in enumerate(order):
        cluster_to_accuracy[cluster] = accuracy[i]
    return cluster_to_accuracy


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

def tf_idf_for_cluster(adata,key):
    df = pd.DataFrame(data=adata.raw.X, index=adata.obs_names, columns=adata.raw.var_names)
    df['cluster'] = adata.obs[key].astype('str').values
    cluster_to_tfidf = {}  # store tfidf score
    cluster_to_exclusive = {}   # store exclusivly expressed genes
    for item in adata.obs[key].cat.categories:
        a = tf_idf_bare_compute(df,item)
        a_names = adata.raw.var_names
        test = pd.Series(data=a, index=a_names)
        test.sort_values(ascending=False, inplace=True)
        # remove artifact genes
        artifact = pd.read_csv('./artifact_genes.txt',sep='\t')
        artifact_genes = set(artifact['genes'].to_list())
        test_pure = test.loc[~test.index.isin(artifact_genes)]
        result = test_pure.iloc[9]
        #print('{0} has {1}'.format(item, result))
        cluster_to_tfidf[item] = result
        cluster_to_exclusive[item] = test[:10].to_dict()
    pd.Series(cluster_to_exclusive,name='genes').to_csv('./scTriangulate_result/exclusive_gene_{}.txt'.format(key),sep='\t')
    return cluster_to_tfidf




def SCCAF_score(adata, key):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix
    # define X and Y and exclude cells whose cluster only have 1 cell
    X = adata.raw.X
    Y = adata.obs[key].values

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
    pd.DataFrame(data=m,index=le.classes_,columns=le.classes_).to_csv('./scTriangulate_result/confusion_SCCAF_{}.txt'.format(key),sep='\t')
    # derive cluster reliability from confusion matrix for each cluster
    numeric2reliable = []  # [0.4,0.5...] length is the number of clusters involved in self-projection
    for i in range(m.shape[0]):
        numeric2reliable.append(m[i, i] / m[i, :].sum())
    cluster_to_SCCAF = {}
    for i in range(len(numeric2reliable)):
        #print('{0} has {1}'.format(le.classes_[i],numeric2reliable[i]))
        cluster_to_SCCAF[le.classes_[i]] = numeric2reliable[i]
    return cluster_to_SCCAF

























# ## code storage
# def tf_idf_all_gene(df,cluster):
#     '''
#     now the df contains all the gene for and an additional column for cluster
#     '''
#     # compute its tf_idf
#     tmp1 = df.loc[df['cluster'] == cluster, :].loc[:,df.columns!='cluster'].values  # (n_cells,n_genes)
#     tf = np.count_nonzero(tmp1,axis=0) / tmp1.shape[0]  # (n_genes,)
#     tf = tf + 1e-5
#     tmp2 = df.loc[:,df.columns!='cluster'].values
#     df_ = np.count_nonzero(tmp2,axis=0) / tmp2.shape[0]  # (n_genes,)
#     df_ = df_ + 1e-5
#     idf = -np.log10(df_)
#     tf_idf_ori = tf * idf  # (n_genes,)
#     # compute background tf_idf
#     ntimes = 1000
#     background = np.empty([ntimes,len(tf_idf_ori)])  # (ntimes,n_genes)
#     for j in range(ntimes):
#         cluster_col = df['cluster'].values
#         np.random.shuffle(cluster_col)
#         df['cluster'] = cluster_col
#         tmp1 = df.loc[df['cluster'] == cluster, :].loc[:, df.columns != 'cluster'].values  # (n_cells,n_genes)
#         tf = np.count_nonzero(tmp1, axis=0) / tmp1.shape[0]  # (n_genes,)
#         tf = tf + 1e-5
#         tmp2 = df.loc[:, df.columns != 'cluster'].values
#         df_ = np.count_nonzero(tmp2, axis=0) / tmp2.shape[0]  # (n_genes,)
#         df_ = df_ + 1e-5
#         idf = -np.log10(df_)
#         tf_idf = tf * idf  # (n_genes,)
#         background[j,:] = tf_idf
#     # compute p_value
#     pval = 1 / ntimes * np.count_nonzero(np.array(background) > tf_idf_ori, axis=0)  # (n_genes,)
#     return tf_idf_ori, pval

# def tf_idf_score(adata,key):
#     big = np.empty([3,len(adata.obs[key].cat.categories),len(adata.raw.var_names)])
#     '''
#     the big ndarray:
#     depth is three layers, ori,pval,flag
#     height is n_cluster
#     width is n_gene
#     '''
#     df = pd.DataFrame(data=adata.raw.X,index=adata.obs_names,columns=adata.raw.var_names)
#     df['cluster'] = adata.obs[key].astype('str').values
#     for i,cluster in enumerate(adata.obs[key].cat.categories):
#         tf_idf_ori, pval = tf_idf_all_gene(df,cluster)
#         flag = pval < 0.05
#         big[0,i,:] = tf_idf_ori
#         big[1,i,:] = pval
#         big[2,i,:] = flag
#     amount = np.count_nonzero(big[2,:,:],axis=1)  # (n_cluster,)
#     print(adata.obs[key].cat.categories)
#     print(amount)