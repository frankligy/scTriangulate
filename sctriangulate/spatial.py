import squidpy as sq
import scanpy as sc
import os,sys
import numpy as np
import pandas as pd
import anndata as ad
import networkx as nx
from .preprocessing import *
from PIL import Image
from json import load


# spatial IO
def read_spatial_data(mode_count='mtx',mode_spatial='visium',mtx_folder=None,spatial_folder=None,spatial_library_id=None,**kwargs):
    '''
    read the spatial data into the memory as adata

    :param mode_count: string, how the spatial count data is present, it can be a folder containing mtx file, or h5 file, or others
    :param mode_spatial: string, how the spatial images and associated files are present, it can be in visium format, or others
    :param mtx_folder: string, if mode_count == 'mtx', specify the folder name
    :param spatial_folder: string, if mode_spatial == 'visium', specifiy the folder name
    :param spatial_library_id, string, when necessary, specify the library_id for the spatial slide
    :param **kwargs: optional keyword arguments will be passed to mtx_to_adata

    Examples::
            id_ = '1160920F'
            adata_spatial = read_spatial_data(mtx_folder='filtered_count_matrices/{}_filtered_count_matrix'.format(id_),
                                              spatial_folder='filtered_count_matrices/{}_filtered_count_matrix/{}_spatial'.format(id_,id_),
                                              spatial_library_id=id_,feature='features.tsv')

    '''
    if mode_count == 'mtx':
        adata_spatial = mtx_to_adata(int_folder=mtx_folder,**kwargs)
    if mode_spatial == 'visium':
        # spatial coordinate
        coords = pd.read_csv(os.path.join(spatial_folder,'tissue_positions_list.csv'), index_col=0, header=None)
        coords.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
        tmp = pd.merge(adata_spatial.obs, coords, how="left", left_index=True, right_index=True)
        tmp.index.name = None
        adata_spatial.obs = tmp
        adata_spatial.obsm['spatial'] = adata_spatial.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values
        adata_spatial.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)
        # uns
        adata_spatial.uns['spatial'] = {}
        adata_spatial.uns['spatial'][spatial_library_id] = {}
        ### uns images
        adata_spatial.uns['spatial'][spatial_library_id]['images'] = {}
        for res in ['hires','lowres']:
            adata_spatial.uns['spatial'][spatial_library_id]['images'][res] = np.array(Image.open(os.path.join(spatial_folder,'tissue_{}_image.png'.format(res))))
        ### uns scalefactors
        f = open(os.path.join(spatial_folder,'scalefactors_json.json'))
        adata_spatial.uns['spatial'][spatial_library_id]['scalefactors'] = load(f)
        ### uns metadata
        adata_spatial.uns['spatial'][spatial_library_id]['metadata'] = {}
    return adata_spatial
    






# implement spatial additional stability metrics cluster level
def cluster_level_spatial_stability(adata,key,method,neighbor_key='spatial_distances',sparse=True,coord_type='generic',n_neighs=6,radius=None,delaunay=False):
    '''
    derive optional stability score in the context of spatial transcriptomics

    :param adata: the Anndata
    :param key: string, the column in obs to derive cluster-level stability score
    :param method: string, which score, support tbe following:

                   # degree_centrality
                   # closeness_centrality
                   # average_clustering 
                   # spread
                   # assortativity
                   # number_connected_components

    :param neighbor_key: string, which obsm key to use for neighbor adjancency, default is spatial_distances
    :param sparse: boolean, whether the adjancency matrix is sparse or not, default is True
    :param coord_type: string, default is generic, passed to sq.gr.spatial_neighbors()
    :param n_neighs: int, default is 6, passed to sq.gr.spatial_neighbors()
    :param radius: float or None (default), passed to sq.gr.spatial_neighbors()
    :param delaunay: boolean, default is False, whether to use delaunay for spatial graph, passed to sq.gr.spatial_neighbors()

    Examples::
        cluster_level_spatial_stability(adata,'cluster',method='centrality')
        cluster_level_spatial_stability(adata,'cluster',method='spread')
        cluster_level_spatial_stability(adata,'cluster',method='assortativity',neighbor_key='spatial_distances',sparse=True)
        cluster_level_spatial_stability(adata,'cluster',method='number_connected_components',neighbor_key='spatial_distances',sparse=True)

    '''
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
        sq.gr.spatial_neighbors(adata,coord_type=coord_type,n_neighs=n_neighs,radius=radius,delaunay=delaunay)
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
        sq.gr.spatial_neighbors(adata,coord_type=coord_type,n_neighs=n_neighs,radius=radius,delaunay=delaunay)
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
def create_spatial_features(adata,mode,coord_type='generic',n_neighs=6,radius=None,delaunay=False,sparse=True,
                            library_id=None,img_key='hires_image',sf_key='tissue_hires_scalef',sd_key='spot_diameter_fullres',feature_types=['summary','texture','histogram'],
                            feature_added_kwargs=[{},{},{}],segmentation_feature=False,segmentation_method='watershed'):
    '''
    Extract spatial features (including spatial coordinates, spatial neighbor graph, spatial image)

    :param adata: the adata to extract features from
    :param mode: string, support:

        # coordinate
        # graph_importance
        # tissue_images

    :param coord_type: string, default is generic, passed to sq.gr.spatial_neighbors()
    :param n_neighs: int, default is 6, passed to sq.gr.spatial_neighbors()
    :param radius: float or None (default), passed to sq.gr.spatial_neighbors()
    :param delaunay: boolean, default is False, whether to use delaunay for spatial graph, passed to sq.gr.spatial_neighbors()
    :param library_id: string or None, when choosing tissue_images, need to retrieve the image from adata.uns
    :param img_key: string, when choosing tissue_images, need to retrieve the image from adata.uns
    :param sf_key: string, when choosing tissue_images, we need to associate image pixel to the spot location, which relies on scalefactor
    :param sd_key: the string, the key of spot_diameter_fullres in adata.uns['spatial'][libarary_id]['scalefactorf']
    :param feature_types: list, default including summary, texture and histogram
    :param feature_added_kwargs: nested list, each element is a dict, containing additional keyword arguments being passed to each feature function in feature_types
    :param segmentation_feature: boolean, whether to extract segmentation feature or not, default is False
    :param segmentation_method: string the segmentation method to use, default is 'watershed'

    '''
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
    elif mode == 'tissue_image':
        # to comply with convention, we assume the img is stored at adata.uns['spatial'][{library_id}]['images'], at certain key slot, dims are (y,x,channels)
        # we assume scalefactor will be adata.uns['spatial'][{library_id}]['scalefactors'], at certain key slot
        # we assume spot_diameter will be adata.uns['spatial'][{library_id}]['scalefactors'], at certain key slot
        # if adding custom function, feature_types = ['summary','texture',histogram','custom'], define one custom function that return all you want
        # then for feature_added_kwargs = [{},{},{},{'func':customized_func,'arg2':4}]
        scalef = adata.uns['spatial'][library_id]['scalefactors'][sf_key]
        img = sq.im.ImageContainer(img=adata.uns['spatial'][library_id]['images'][img_key],layer='image',dims=('y','x','channels'),scale=scalef)
        adata_copy = adata.copy()
        features_df_list = []
        for feature,added_kwargs in zip(feature_types,feature_added_kwargs):
            sq.im.calculate_image_features(adata_copy, img, features=feature, features_kwargs={feature: added_kwargs},key_added="{}_features".format(feature), show_progress_bar=True)
            features_df_list.append(adata_copy.obsm["{}_features".format(feature)])
        if segmentation_feature:
            sq.im.segment(img=img, layer="image", layer_added="segmented", method=segmentation_method)
            sq.im.calculate_image_features(adata_copy,img,layer="image",features="segmentation",key_added="segmentation_features",features_kwargs={"segmentation": {"label_layer": "segmented"}},mask_circle=True)
            features_df_list.append(adata_copy.obsm['segmentation_features'])
        spatial_adata = ad.AnnData(pd.concat(features_df_list,axis=1))
    spatial_adata.obsm['spatial'] = adata.obsm['spatial']
    return spatial_adata



