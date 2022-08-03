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
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

# for publication ready figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# spatial IO
def read_spatial_data(mode_count='mtx',mode_spatial='visium',mtx_folder=None,spatial_folder=None,spatial_library_id=None,**kwargs):
    '''
    read the spatial data into the memory as adata
    
    :param mode_count: string, how the spatial count data is present, it can be a folder containing mtx file, or h5 file, or others
    :param mode_spatial: string, how the spatial images and associated files are present, it can be in visium format, or others
    :param mtx_folder: string, if mode_count == 'mtx', specify the folder name
    :param spatial_folder: string, if mode_spatial == 'visium', specifiy the folder name
    :param spatial_library_id: string, when necessary, specify the library_id for the spatial slide
    :param kwargs: optional keyword arguments will be passed to mtx_to_adata

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
    


# visualize deconvolution
def plot_one_dot(ax,prop,x,y,scale_factor,circle_diameter,colors):
    previous = 0
    for i, (p, c) in enumerate(zip(prop, colors)):
        this = 2 * np.pi * p + previous
        x_  = [0] + np.cos(np.linspace(previous, this, 50)).tolist() + [0]
        y_  = [0] + np.sin(np.linspace(previous, this, 50)).tolist() + [0]
        xy_ = np.column_stack([x_, y_])
        s = np.abs(xy_).max()
        previous = this
        ax.scatter([x * scale_factor], [y * scale_factor], marker=xy_, s=s**2 * circle_diameter, facecolor=c,lw=0)
    return ax

def plot_deconvolution(adata,decon,fraction=0.1,size=0.5,library_id=None,alpha=0.5,outdir='.'):
    '''
    Visualize the deconvolution result as pie chart, serving as a complement function on top of scanpy.pl.spatial where each dot now is a pie chart

    :param adata: AnnData, requires img, scale_factor, spot_diameter in the canonical slot based on squidpy convention, you can get it either using squidpy.read_visium or sctriangulate.spatial.read_spatial_data
    :param decon: A dataframe whose index is spot barcode, column is the cell types, value represents cell type proportion (each row sums to 1)
    :param fraction: float, plot that fraction of spot on the image, default is 0.1
    :param size: float, default is 0.5, adjust it based on visual effect
    :param library_id: None or string, if None, automatically extract from adata.uns['spatial'], or particularly specify
    :param alpha: float, default is 0.5, the transparency of the underlying image
    :param outdir: string, the output dir for saved image

    Example::

        decon = pd.read_csv('inputs/decon_results_A/cell2location_prop_processed.txt',index_col=0,sep='\t')
        plot_deconvolution(adata_spatial,decon,fraction=0.5,alpha=0.3)

    .. image:: ./_static/decon.png
        :height: 400px
        :width: 550px
        :align: center
        :target: target  

    '''
    decon = decon.loc[adata.obs_names,:]
    adata.obsm['decon'] = decon
    adata = adata.copy()
    sc.pp.subsample(adata,fraction=fraction)
    if library_id is None:
        library_id = list(adata.uns['spatial'].keys())[0]
        spot_size = adata.uns['spatial'][library_id]['scalefactors']['spot_diameter_fullres']
        scale_factor = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
        img = adata.uns['spatial'][library_id]['images']['hires']
        circle_diameter = spot_size * scale_factor * size
        colors_set = colors_for_set(adata.obsm['decon'].columns.tolist())
        colors_list = list(colors_set.values())
        fig, ax = plt.subplots()
        for prop, (x, y) in tqdm(zip(adata.obsm['decon'].values, adata.obsm['spatial']), total=adata.shape[0]):
            ax = plot_one_dot(ax, prop, x, y, scale_factor, circle_diameter, colors_list)
        ax.imshow(X=img, aspect='equal', alpha=alpha, origin='upper')
        ax.set_xticks([])
        ax.set_yticks([])
        import matplotlib.patches as mpatches
        ax.legend(handles=[mpatches.Patch(color=i) for i in colors_list], labels=list(colors_set.keys()),
                  loc='upper left', bbox_to_anchor=[1, 1], frameon=False)
        plt.savefig(os.path.join(outdir,'deconvolution_fraction_{}_size_{}_alpha_{}.pdf'.format(fraction,size,alpha)), bbox_inches='tight')
        plt.close()






# implement spatial additional stability metrics cluster level
def cluster_level_spatial_stability(adata,key,method,neighbor_key='spatial_distances',sparse=True,coord_type='generic',n_neighs=6,radius=None,delaunay=False):
    '''
    derive optional stability score in the context of spatial transcriptomics

    :param adata: the Anndata
    :param key: string, the column in obs to derive cluster-level stability score
    :param method: string, which score, support tbe following:

                * degree_centrality
                * closeness_centrality
                * average_clustering 
                * spread
                * assortativity
                * number_connected_components

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

        * **coordinate**: feature derived from pure spatial coordinates
        * **graph_importance**: feature derived from the spatial neighbor graph 
        * **tissue_images**: feature derived from assciated tissue images (H&E, fluorescent)

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


def multicellular_spatial_structure(adata_spatial,coord_type='grid',n_neighrs=6,n_rings=1,spatial_community_resolution=8,spatial_community_min_cells=10,
                                    plot=True,outdir='.',mode='spot_cell_type_proportion'):
    '''
    The most important analysis is to define **spatial program**, meaning spatial region that are proximal and somehow transcriptomically similar. This function
    provide a very flexible way to define `spatial program` by first derive spatial community solely based on spatial coordinates, then we merge spatial community into
    spatial cluster/structure such that spatial communities sharing similar gene expression or cell phenotype profile will be clustered tegether, the resultant spatial cluster
    represents the biologically meaningful multicellular spatial structure, it is inspired by Figure 4 in `this paper <https://www.nature.com/articles/s41588-022-01041-y>`_.

    :param adata_spatial: the Anndata with spatial coordinate
    :param coord_type: either grid or generic, passed to sq.gr.spatial_neighbors function
    :param n_neighrs: int, default is 6, passed to sq.gr.spatial_neighbors function
    :param n_rings: int, default is 1, passed to sq.gr.spatial_neighbors function
    :param spatial_community_resolution: float, default is 8, passed to nx.algorithms.community.greedy_modularity_communities, high value means smaller clusters
    :param spatial_community_min_cells: int, the minimum number of cells/spots that a spatial community needs to have
    :param plot: boolean, whether to plot spatial community and spatial cluster
    :param outdir: string, the path to the output dir
    :param mode: string, which transcriptomical profile to consider, it can be 'spot_cell_type_proportion' coming out of spatial deconvolution methods or spot gene expression

    Examples::

        multicellular_spatial_structure(adata_spatial,outdir='result')

    '''
    adata_spatial = adata.spatial.copy()
    # step1: build spatial community
    sq.gr.spatial_neighbors(adata_spatial,coord_type=coord_type,n_neighs=n_neighbors,n_rings=n_rings)
    G = nx.from_scipy_sparse_matrix(adata_spatial.obsp['spatial_distances'])
    index2barcode = pd.Series(index=list(G.nodes()),data=adata_spatial.obs_names.tolist()).to_dict()
    partitions = nx.algorithms.community.greedy_modularity_communities(G,resolution=spatial_community_resolution)
    mapping = {}
    for i,par in enumerate(partitions):
        if len(par) >= spatial_community_min_cells:
            for node in par:
                mapping[index2barcode[node]] = i
    adata_spatial.obs['community'] = adata_spatial.obs_names.map(mapping).fillna('unknown').values
    adata_spatial = adata_spatial[adata_spatial.obs['community']!='unknown',:].copy()
    adata_spatial.obs['community'] = adata_spatial.obs['community'].astype('int').astype('category')
    if plot:
        sc.pl.spatial(adata_spatial,color='community')
        plt.savefig(os.path.join(outdir,'spatial_community_plot.pdf'))
        plt.close()
    # step2: integrate other information and cluster communities
    if mode == 'spot_cell_type_proportion':
        # assuming the cell type proportion for each spot is in adata.X
        decon = adata.spatial.to_df()
        interactions = list(product(decon.columns.tolist(),decon.columns.tolist()))
        data_list = []
        for i,par in tqdm(enumerate(partitions),total=len(partitions)):
            if len(par) >= spatial_community_min_cells:
                G_sub = G.subgraph(nodes=par).copy()
                nodes = [index2barcode[node] for node in G_sub.nodes()]
                feature = []
                for inter in interactions:
                    v = 0
                    tmp = adata_decon[nodes,inter]
                    for edge in G_sub.edges():
                        c1 = index2barcode[edge[0]]
                        c2 = index2barcode[edge[1]]
                        v += tmp[[c1,c2],:].X.mean()
                    v /= len(G_sub.edges())
                    feature.append(v)
                data_list.append(pd.Series(index=interactions,data=feature,name=i))
        df_subgraph = pd.concat(data_list,axis=1).T.astype('float')
        model = AgglomerativeClustering(n_clusters=10)
        model.fit(df_subgraph.values)
        mi_r = pd.MultiIndex.from_arrays([['subgraph_{}'.format(i) for i in df_subgraph.index],['cluster_{}'.format(i) for i in model.labels_]])
        mi_c = pd.MultiIndex.from_tuples(list(df_subgraph.columns))
        df_subgraph.index = mi_r
        df_subgraph.columns = mi_c
        df_subgraph.to_csv(os.path.join(outdir,'df_subgraph.txt'),sep='\t')
        if plot:
            mapping2 = pd.Series(index=df_subgraph.index,data=model.labels_).to_dict()
            adata_spatial.obs['cluster'] = adata_spatial.obs['community'].map(mapping2).values
            adata_spatial.obs['cluster'] = adata_spatial.obs['cluster'].astype('category')
            sc.pl.spatial(adata_decon,color='cluster')
            plt.savefig(os.path.join(outdir,'spatial_structure_plot.pdf'))
            plt.close()            

        return df_subgraph

