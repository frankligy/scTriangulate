import squidpy as sq
import scanpy as sc
import os,sys
import numpy as np
import pandas as pd
import anndata as ad
import networkx as nx



# implement spatial additional stability metrics cluster level
def cluster_level_spatial_stability(adata,key,method,neighbor_key=None,sparse=None,coord_type='generic',n_neighs=6,radius=None,delaunay=False):
    if method == 'degree_centrality' or method == 'closeness_centrality' or method == 'average_clustering':
        sq.gr.spatial_neighbors(adata,coord_type=coord_type,n_neighs=n_neighs,radius=radius,delaunay=delaunay)
        sq.gr.centrality_scores(adata, cluster_key=key)
        df = adata.uns['{}_centrality_scores'.format(key)]
        mapping = df[method].to_dict()
    elif method =='spread':
        sq.gr.ripley(adata, cluster_key=key, mode='L')
        df = adata.uns['{}_ripley_L'.format(key)]['L_stat']
        q = 0.5
        mapping = {}
        for c, subdf in df.groupby(by=key):
            t = subdf.shape[0]
            s = round(t * q)
            v = subdf.iloc[s, :]['stats']
            mapping[c] = v
    elif method == 'assortativity':
        adjacency_matrix = adata.obsp[neighbor_key]
        if sparse:
            G = nx.from_scipy_sparse_matrix(adjacency_matrix)
        else:
            G = nx.from_numpy_matrix(adjacency_matrix)
        index2cluster = pd.Series(index=np.arange(adata.shape[0]), data=adata.obs[key].values).to_dict()
        all_cluster = adata.obs[key].cat.categories.tolist()
        all_index = np.arange(len(all_cluster))
        cluster2order = pd.Series(index=all_cluster, data=all_index).to_dict()
        nx.set_node_attributes(G, index2cluster, 'cluster')
        mix_mat = nx.attribute_mixing_matrix(G, attribute='cluster', mapping=cluster2order)
        mapping = {}
        for c, o in cluster2order.items():
            self_mix = mix_mat[o, o] / mix_mat[o, :].sum()
            mapping[c] = o
    elif method == 'number_connected_components':
        adjacency_matrix = adata.obsp[neighbor_key]
        if sparse:
            G = nx.from_scipy_sparse_matrix(adjacency_matrix)
        else:
            G = nx.from_numpy_matrix(adjacency_matrix)
        adata.obs['index'] = np.arange(adata.shape[0])
        cluster2index = adata.obs.loc[:, [key, 'index']].set_index(keys='index').reset_index().groupby(by=key)['index'].apply(lambda x: x.tolist()).to_dict()
        mapping = {}
        for c, indices in cluster2index.items():
            subgraph = G.subgraph(nodes=indices)
            n = nx.number_connected_components(subgraph)
            mapping[c] = n
    return mapping


# implement spatial cell-level metrics
def create_spatial_features(adata,mode,coord_type='generic',n_neighs=6,radius=None,delaunay=False,sparse=True):
    if mode == 'coordinate':
        spatial_adata = ad.AnnData(X=adata.obsm['spatial'],var=pd.DataFrame(index=['spatial_x','spatial_y']),obs=pd.DataFrame(index=adata.obs_names))
    elif mode == 'graph_importance':
        adata_copy = adata.copy()
        sq.gr.spatial_neighbors(adata_copy,coord_type=coord_type,n_neighs=n_neighs,radius=radius,delaunay=delaunay)
        adjacency_matrix = adata_copy.obsp['spatial_distances']
        if sparse:
            G = nx.from_scipy_sparse_matrix(adjacency_matrix)
        else:
            G = nx.from_numpy_matrix(adjacency_matrix)
        values_dict = {}
        degree_dict = dict(G.degree(weight='weight')); values_dict['degree'] = degree_dict
        cc_dict = nx.clustering(G, weight='weight'); values_dict['clustering_coeffient'] = cc_dict
        pr_dict = nx.pagerank(G, weight='weight'); values_dict['pagerank_score'] = pr_dict
        for k,v in values_dict.items():
            adata_copy.obs[k] = list(v.values())
        spatial_adata = ad.AnnData(adata_copy.obs[list(values_dict.keys())])
    spatial_adata.obsm['spatial'] = adata.obsm['spatial']
    return spatial_adata



