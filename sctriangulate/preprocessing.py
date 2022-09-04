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
import scanpy as sc
import anndata as ad
from scipy.io import mmread,mmwrite
from scipy.sparse import csr_matrix,issparse
import matplotlib as mpl
from functools import reduce
from sklearn.decomposition import PCA
import umap

from sctriangulate.colors import *



# for publication ready figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def sctriangulate_preprocessing_setting(backend='Agg',png=False):
    # change the backend
    mpl.use(backend)
    if png:
        # for publication and super large dataset
        mpl.rcParams['savefig.dpi'] = 600
        mpl.rcParams['figure.dpi'] = 600


def small_txt_to_adata(int_file,gene_is_index=True):
    '''
    given a small dense expression (<2GB) txt file, load them into memory as adata, and also make sure the X is sparse matrix.

    :param int_file: string, path to the input txt file, delimited by tab
    :param gene_is_index: boolean, whether the gene/features are the index.

    :return: AnnData

    Exmaples::

        from sctriangulate.preprocessing import small_txt_to_adata
        adata= = small_txt_to_adata('./input.txt',gene_is_index=True)

    '''
    df = pd.read_csv(int_file,sep='\t',index_col=0)
    if gene_is_index:
        adata = ad.AnnData(X=csr_matrix(df.values.T),var=pd.DataFrame(index=df.index.values),obs=pd.DataFrame(index=df.columns.values))
    else:
        adata = ad.AnnData(X=csr_matrix(df.values),var=pd.DataFrame(index=df.columns.values),obs=pd.DataFrame(index=df.index.values))
    adata.var_names_make_unique()
    adata.X = csr_matrix(adata.X)
    return adata


def large_txt_to_mtx(int_file,out_folder,gene_is_index=True,type_convert_to='int16'):  # whether the txt if gene * cell
    '''
    Given a large txt dense expression file, convert them to mtx file on cluster to facilitate future I/O

    :param int_file: string, path to the intput txt file, delimited by tab
    :param out_folder: string, path to the output folder where the mtx file will be stored
    :param gene_is_index: boolean, whether the gene/features is the index in the int_file.
    :param type_convert_to: string, since it is a large dataframe, need to read in chunk, to accelarate it and reduce the memory footprint,
                            we convert it to either 'int16' if original data is count, or 'float32' if original data is normalized data.
    
    Examples::

        from sctriangulate.preprocessing import large_txt_to_mtx
        large_txt_to_mtx(int_file='input.txt',out_folder='./data',gene_is_index=False,type_convert_to='float32')
    ''' 
    reader = pd.read_csv(int_file,sep='\t',index_col=0,chunksize=1000)
    store = []
    for chunk in reader:
        tmp = chunk.astype(type_convert_to)
        store.append(tmp)
    data = pd.concat(store)
    print(data.shape)
    '''save as mtx, now!!!'''
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    if gene_is_index:
        data.index.to_series().to_csv(os.path.join(out_folder,'genes.tsv'),sep='\t',header=None,index=None)
        data.columns.to_series().to_csv(os.path.join(out_folder,'barcodes.tsv'),sep='\t',header=None,index=None)
        mmwrite(os.path.join(out_folder,'matrix.mtx'),csr_matrix(data.values))
    else:
        data.columns.to_series().to_csv(os.path.join(out_folder,'genes.tsv'),sep='\t',header=None,index=None)
        data.index.to_series().to_csv(os.path.join(out_folder,'barcodes.tsv'),sep='\t',header=None,index=None)
        mmwrite(os.path.join(out_folder,'matrix.mtx'),csr_matrix(data.values.T))        


def mtx_to_adata(int_folder,gene_is_index=True,feature='genes.tsv',feature_col='index',barcode='barcodes.tsv',barcode_col='index',matrix='matrix.mtx'):  # whether the mtx file is gene * cell
    '''
    convert mtx file to adata in RAM, make sure the X is sparse.

    :param int_folder: string, folder where the mtx files are stored.
    :param gene_is_index: boolean, whether the gene is index.
    :param feature: string, the name of the feature file, default for rna is genes.tsv
    :param feature_col: 'index' as index, or a int (which column, python is zero based) to use in your feature.tsv as feature
    :param barcode: string, the name of the barcode file, default is barcodes.tsv
    :param barcode_col: 'index' as index, or a int (which column, python is zero based) to use in your barcodes.tsv as barcode
    :param matrix: string, the name of the sparse matrix

    :return: AnnData

    Examples::

        from sctriangulate.preprocessing import mtx_to_adata
        mtx_to_adata(int_folder='./data',gene_is_index=False,feature='genes')

    '''
    if feature_col == 'index':
        gene = pd.read_csv(os.path.join(int_folder,feature),sep='\t',index_col=0,header=None).index
    else:
        gene = pd.read_csv(os.path.join(int_folder,feature),sep='\t',index_col=0,header=None)[feature_col]
    if barcode_col == 'index':
        cell = pd.read_csv(os.path.join(int_folder,barcode),sep='\t',index_col=0,header=None).index
    else:
        cell = pd.read_csv(os.path.join(int_folder,barcode),sep='\t',index_col=0,header=None)[barcode_col]
    value = csr_matrix(mmread(os.path.join(int_folder,matrix)))
    if gene_is_index:
        value = value.T
        adata = ad.AnnData(X=value,obs=pd.DataFrame(index=cell),var=pd.DataFrame(index=gene))
    else:
        adata = ad.AnnData(X=value,obs=pd.DataFrame(index=cell),var=pd.DataFrame(index=gene))
    adata.var.index.name = None
    adata.var_names_make_unique()
    return adata


def mtx_to_large_txt(int_folder,out_file,gene_is_index=False):
    '''
    convert mtx back to large dense txt expression dataframe.

    :param int_folder: string, path to the input mtx folder.
    :param out_file: string, path to the output txt file.
    :param gene_is_index: boolean, whether the gene is the index.

    Examples::

        from sctriangulate.preprocessing import mtx_to_large_txt
        mtx_to_large_txt(int_folder='./data',out_file='input.txt',gene_is_index=False)

    '''
    gene = pd.read_csv(os.path.join(int_folder,'genes.tsv'),sep='\t',index_col=0,header=None).index
    cell = pd.read_csv(os.path.join(int_folder,'barcodes.tsv'),sep='\t',index_col=0,header=None).index
    value = mmread(os.path.join(int_folder,'matrix.mtx')).toarray()
    if gene_is_index:
        data = pd.DataFrame(data=value,index=gene,columns=cell)
    else:
        data = pd.DataFrame(data=value.T,index=cell,columns=gene)
    data.to_csv(out_file,sep='\t',chunksize=1000)


def adata_to_mtx(adata,gene_is_index=True,var_column=None,obs_column=None,outdir='data'):
    '''
    write a adata to mtx file

    :param adata: AnnData, the adata to convert
    :param gene_is_index: boolean, for the resultant mtx, will gene be the index or feature is the index
    :param var_column: list, the var columns to write the genes.tsv, None will write all available columns
    :param obs_column: list, the obs columns to write the barcodes.tsv, None will write all available columns
    :param outdir: string, the name of the mtx folder, default is data

    Example::

        adata_to_mtx(adata,True,None,None,'data')
    '''
    # create folder if not exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # write genes.tsv
    if var_column is None:
        var = adata.var_names.to_series()
    else:
        var = adata.var[var_column]
    var.to_csv(os.path.join(outdir,'genes.tsv'),sep='\t',header=None,index=None)
    # write barcodes.tsv
    if obs_column is None:
        obs = adata.obs_names.to_series()
    else:
        obs = adata.obs[obs_column]
    obs.to_csv(os.path.join(outdir,'barcodes.tsv'),sep='\t',header=None,index=None)
    # write matrix.mtx
    if not gene_is_index:
        mmwrite(os.path.join(outdir,'matrix.mtx'),make_sure_mat_sparse(adata.X))
    else:
        mmwrite(os.path.join(outdir,'matrix.mtx'),make_sure_mat_sparse(adata.X).transpose())



def add_azimuth(adata,result,name='predicted.celltype.l2'):
    '''
    a convenient function if you have azimuth predicted labels in hand, and want to add the label to the adata.

    :param adata: AnnData
    :param result: string, the path to the 'azimuth_predict.tsv' file
    :param name: string, the column name where the user want to transfer to the adata.

    Examples::

        from sctriangulate.preprocessing import add_azimuth
        add_azimuth(adata,result='./azimuth_predict.tsv',name='predicted.celltype.l2')

    '''
    azimuth = pd.read_csv(result,sep='\t',index_col=0)
    azimuth_map = azimuth[name].to_dict()
    azimuth_prediction = azimuth['{}.score'.format(name)].to_dict()
    azimuth_mapping = azimuth['mapping.score'].to_dict()
    adata.obs['azimuth'] = adata.obs_names.map(azimuth_map).values
    adata.obs['prediction_score'] = adata.obs_names.map(azimuth_prediction).values
    adata.obs['mapping_score'] = adata.obs_names.map(azimuth_mapping).values

def add_annotations(adata,inputs,cols_input,index_col=0,cols_output=None,sep='\t',kind='disk'):
    '''
    Adding annotations from external sources to the adata

    :param adata: Anndata
    :param inputs: string, path to the txt file where the barcode to cluster label information is stored.
    :param cols_input: list, what columns the users want to transfer to the adata.
    :param index_col: int, for the input, which column will serve as the index column
    :param cols_output: list, corresponding to the cols_input, how these columns will be named in the adata.obs columns
    :param sep: string, default is tab
    :param kind: a string, either 'disk', or 'memory', disk means the input is the path to the text file, 'memory' means the input is the
                variable name in the RAM that represents the dataframe

    Examples::

        from sctriangulate.preprocessing import add_annotations
        add_annotations(adata,inputs='./annotation.txt',cols_input=['col1','col2'],index_col=0,cols_output=['annotation1','annontation2'],kind='disk')
        add_annotations(adata,inputs=df,cols_input=['col1','col2'],index_col=0,cols_output=['annotation1','annontation2'],kind='memory')

    '''
    # means a single file such that one column is barcodes, annotations are within other columns
    if kind == 'disk':
        annotations = pd.read_csv(inputs,sep=sep,index_col=index_col).loc[:,cols_input]
    elif kind == 'memory':   # index_col will be ignored
        annotations = inputs.loc[:,cols_input]
    mappings = []
    for col in cols_input:
        mapping = annotations[col].to_dict()
        mappings.append(mapping)
    if cols_output is None:
        for i,col in enumerate(cols_input):
            adata.obs[col] = adata.obs_names.map(mappings[i]).fillna('Unknown').values
            adata.obs[col] = adata.obs[col].astype('str').astype('category')
    else:
        for i in range(len(cols_input)):
            adata.obs[cols_output[i]] = adata.obs_names.map(mappings[i]).fillna('Unknown').values
            adata.obs[cols_output[i]] = adata.obs[cols_output[i]].astype('str').astype('category')


def add_umap(adata,inputs,mode,cols=None,index_col=0,key='X_umap'):
    '''
    if umap embedding is pre-computed, add it back to adata object.

    :param adata: Anndata
    :param inputs: string, path to the the txt file where the umap embedding was stored.
    :param mode: string, valid value 'pandas_disk', 'pandas_memory', 'numpy'

        * **pandas_disk**: the `inputs` argument should be the path to the txt file
        * **pandas_memory**: the `inputs` argument should be the name of the pandas dataframe in the program, inputs=df
        * **numpy**, the `inputs` argument should be a 2D ndarray contains pre-sorted (same order as barcodes in adata) umap coordinates

    :param cols: list, what columns contain umap embeddings
    :param index_col: int, which column will serve as the index column.
    :param key: string, the key name to add to obsm

    Examples::

        from sctriangulate.preprocessing import add_umap
        add_umap(adata,inputs='umap.txt',mode='pandas_disk',cols=['umap1','umap2'],index_col=0)

    '''
    # make sure cols are [umap_x, umap_y]
    if mode == 'pandas_disk':
        df = pd.read_csv(inputs,sep='\t',index_col=index_col)
        umap_x = df[cols[0]].to_dict()
        umap_y = df[cols[1]].to_dict()
        adata.obs['umap_x'] = adata.obs_names.map(umap_x).values
        adata.obs['umap_y'] = adata.obs_names.map(umap_y).values
        adata.obsm[key] = adata.obs.loc[:,['umap_x','umap_y']].values
        adata.obs.drop(columns=['umap_x','umap_y'],inplace=True)
    elif mode == 'pandas_memory':
        df = inputs
        umap_x = df[cols[0]].to_dict()
        umap_y = df[cols[1]].to_dict()
        adata.obs['umap_x'] = adata.obs_names.map(umap_x).values
        adata.obs['umap_y'] = adata.obs_names.map(umap_y).values
        adata.obsm[key] = adata.obs.loc[:,['umap_x','umap_y']].values
        adata.obs.drop(columns=['umap_x','umap_y'],inplace=True)
    elif mode == 'numpy':  # assume the order is correct
        adata.obsm[key] = inputs

def doublet_predict(adata):  # gave RNA count or log matrix
    '''
    wrapper function for running srublet, a new column named 'doublet_scores' will be added to the adata

    :param adata: Anndata

    :return: dict

    Examples::

        from sctriangulate.preprocessing import doublet_predict
        mapping = doublet_predict(old_adata)

    '''
    from scipy.sparse import issparse
    import scrublet as scr
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    counts_matrix = adata.X
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, min_cells=1)
    adata.obs['doublet_scores'] = doublet_scores
    return adata.obs['doublet_scores'].to_dict()

def make_sure_adata_writable(adata,delete=False):
    '''
    maks sure the adata is able to write to disk, since h5 file is stricted typed, so no mixed dtype is allowd.
    this function basically is to detect the column of obs/var that are of mixed types, and delete them.

    :param adata: Anndata
    :param delete: boolean, False will just print out what columns are mixed type, True will automatically delete those columns

    :return: Anndata

    Examples::

        from sctriangulate.preprocessing import make_sure_adata_writable
        make_sure_adata_writable(adata,delete=True)


    '''
    # check index, can not have name
    var_names = adata.var_names
    obs_names = adata.obs_names
    var_names.name = None
    obs_names.name = None
    adata.var_names = var_names
    adata.obs_names = obs_names

    # make sure column is pure type, basically, if mixed tyep, delete the column, and print out the delete columns
    # go to: https://github.com/theislab/scanpy/issues/1866

    var = adata.var
    obs = adata.obs
    
    for col in var.columns:
        if var[col].dtypes == 'O':
            all_type = np.array([type(item) for item in var[col]])
            first = all_type[0]
            if (first==all_type).all() and first == str:   # object, but every item is str
                continue
            else:   # mixed type
                print('column {} in var will be deleted, because mixed types'.format(col))
                if delete:
                    adata.var.drop(columns=[col],inplace=True)
    
    for col in obs.columns:
        if obs[col].dtypes == 'O':
            all_type = np.array([type(item) for item in obs[col]])
            first = all_type[0]
            if (first==all_type).all() and first == str:   # object, but every item is str
                continue
            else:   # mixed type
                print('column {} in obs will be deleted, because mixed types'.format(col))
                if delete:
                    adata.obs.drop(columns=[col],inplace=True)

    return adata


def scanpy_recipe(adata,species='human',is_log=False,resolutions=[1,2,3],modality='rna',umap=True,save=True,pca_n_comps=None,n_top_genes=3000):
    '''
    Main preprocessing function. Run Scanpy normal pipeline to achieve Leiden clustering with various resolutions across multiple modalities.

    :param adata: Anndata
    :param species: string, 'human' or 'mouse'
    :param is_log: boolean, whether the adata.X is count or normalized data.
    :param resolutions: list, what leiden resolutions the users want to obtain.
    :param modality: string, valid values: 'rna','adt','atac', 'binary'[mutation data, TCR data, etc]
    :param umap: boolean, whether to compute umap embedding.
    :param save: boolean, whether to save the obtained adata object with cluster label information in it.
    :param pca_n_comps: int, how many PCs to keep when running PCA. Suggestion: RNA (30-50), ADT (15), ATAC (100)
    :param n_top_genes: int, how many features to keep when selecting highly_variable_genes. Suggestion: RNA (3000), ADT (ignored), ATAC (50000-100000)

    :return: Anndata

    Examples::

        from sctriangulate.preprocessing import scanpy_recipe
        # rna
        adata = scanpy_recipe(adata,is_log=False,resolutions=[1,2,3],modality='rna',pca_n_comps=50,n_top_genes=3000)
        # adt
        adata = scanpy_recipe(adata,is_log=False,resolutions=[1,2,3],modality='adt',pca_n_comps=15)
        # atac
        adata = scanpy_recipe(adata,is_log=False,resolutions=[1,2,3],modality='atac',pca_n_comps=100,n_top_genes=100000)
        # binary
        adata = scanpy_recipe(adata,resolutions=[1,2,3],modality='binary')

    '''
    adata.var_names_make_unique()
    # normal analysis
    if modality == 'rna':
        if not is_log:   # count data
            if species == 'human':
                adata.var['mt'] = adata.var_names.str.startswith('MT-')
            elif species == 'mouse':
                adata.var['mt'] = adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
            sc.pp.normalize_total(adata,target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=n_top_genes)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
            sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            # put raw back to X, and make sure it is sparse matrix
            adata = adata.raw.to_adata()
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))
 

        else:   # log(1+x) and depth normalized data
            if species == 'human':
                adata.var['mt'] = adata.var_names.str.startswith('MT-')
            elif species == 'mouse':
                adata.var['mt'] = adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=n_top_genes)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
            sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            # put raw back to X, and make sure it is sparse matrix
            adata = adata.raw.to_adata()
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))

    elif modality == 'atac':
        if not is_log:
            sc.pp.calculate_qc_metrics(adata,percent_top=None,inplace=True,log1p=False)
            sc.pp.normalize_total(adata,target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=n_top_genes)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            #sc.pp.scale(adata,max_value=10)   # because in episcanpy toturial, it seems to be ignored
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            adata = adata.raw.to_adata()
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))
        else:
            sc.pp.calculate_qc_metrics(adata,percent_top=None,inplace=True,log1p=False)
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=n_top_genes)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            #sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            adata = adata.raw.to_adata()
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))

    elif modality == 'adt':
        if not is_log:
            sc.pp.calculate_qc_metrics(adata,percent_top=None,inplace=True,log1p=False)
            adata.X = make_sure_mat_sparse(Normalization.CLR_normalization(make_sure_mat_dense(adata.X)))
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))
        else:
            sc.tl.pca(adata,n_comps=pca_n_comps)
            sc.pp.neighbors(adata)
            for resolution in resolutions:
                sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
            if umap:
                sc.tl.umap(adata)
            if not issparse(adata.X):
                adata.X = csr_matrix(adata.X)
            if save:
                resolutions = '_'.join([str(item) for item in resolutions])
                adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))

    elif modality == 'binary':  # mutation 
        #sc.tl.pca(adata,n_comps=pca_n_comps)
        sc.pp.neighbors(adata,use_rep='X',metric='jaccard')
        for resolution in resolutions:
            sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
        if umap:
            sc.tl.umap(adata)
        if not issparse(adata.X):
            adata.X = csr_matrix(adata.X)
        if save:
            resolutions = '_'.join([str(item) for item in resolutions])
            adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,umap))

    elif modality == 'spatial':
        sc.pp.scale(adata)
        sc.pp.neighbors(adata)
        for resolution in resolutions:
            sc.tl.leiden(adata,resolution=resolution,key_added='sctri_{}_leiden_{}'.format(modality,resolution))
        if save:
            resolutions = '_'.join([str(item) for item in resolutions])
            adata.write('adata_after_scanpy_recipe_{}_{}_umap_{}.h5ad'.format(modality,resolutions,False))        

    return adata
    




def concat_rna_and_other(adata_rna,adata_other,umap,umap_key,name,prefix):
    '''
    concatenate rna adata and another modality's adata object

    :param adata_rna: AnnData
    :param adata_other: Anndata
    :param umap: string, whose umap to use, either 'rna' or 'other'
    :param name: string, the name of other modality, for example, 'adt' or 'atac'
    :param prefix: string, the prefix added in front of features from other modality, by scTriangulate convertion, adt will be 'AB_', atac will be ''.

    :return adata_combine: Anndata

    Examples::

        from sctriangulate.preprocessing import concat_rna_and_other
        concat_rna_and_other(adata_rna,adata_adt,umap='rna',name='adt',prefix='AB_') 

    '''
    adata_rna = adata_rna.copy()
    adata_other = adata_other.copy()
    # remove layers, [!obsm], varm, obsp, varp, raw
    for adata in [adata_rna,adata_other]:
        del adata.layers
        del adata.varm
        del adata.obsp
        del adata.varp
        del adata.raw
    adata_other = adata_other[adata_rna.obs_names,:]  # make sure the obs order is the same
    adata_other.var_names = [prefix + item for item in adata_other.var_names]
    adata_combine = ad.concat([adata_rna,adata_other],axis=1,join='outer',merge='first',label='modality',keys=['rna','{}'.format(name)])
    if umap == 'rna':
        adata_combine.obsm[umap_key] = adata_rna.obsm[umap_key]
    elif umap == 'other':
        adata_combine.obsm[umap_key] = adata_other.obsm[umap_key]
    if not issparse(adata_combine.X):
        adata_combine.X = csr_matrix(adata_combine.X)
    return adata_combine


def nca_embedding(adata,nca_n_components,label,method,max_iter=50,plot=True,save=True,format='pdf',legend_loc='on data',n_top_genes=None,hv_features=None,add_features=None):
    '''
    Doing Neighborhood component ananlysis (NCA), so it is a supervised PCA that takes the label from the annotation, and try to generate a UMAP 
    embedding that perfectly separate the labelled clusters.

    :param adata: the Anndata
    :param nca_n_components: recommend to be 10 based on `Ref <https://www.nature.com/articles/s41586-021-03969-3>`_
    :param label: string, the column name which contains the label information
    :param method: either 'umap' or 'tsne'
    :param max_iter: for the NCA, default is 50, it is generally good enough
    :param plot: whether to plot the umap/tsne or not
    :param save: whether to save the plot or not
    :param format: the saved format, default is 'pdf'
    :param legend_loc: 'on data' or 'right margin'
    :param n_top_genes: how many hypervariable genes to choose for NCA, recommended 3000 or 5000, default is None, means there will be other features to add, multimodal setting
    :param hv_features: a list contains the user-supplied hypervariable genes/features, in multimodal setting, this can be [rna genes] + [ADT protein]
    :param add_features: this should be another adata contains features from other modalities, or None means just for RNA

    Example::

        from sctriangulate.preprocessing import nca_embedding
        # only RNA
        nca_embedding(adata,nca_n_components=10,label='annotation1',method='umap',n_top_genes=3000)
        # RNA + ADT
        # list1 contains [gene features that are variable] and [ADT features that are variable]
        nca_embedding(adata_rna,nca_n_components=10,label='annotation1',method='umap',n_top_genes=3000,hv_features=list1, add_features=adata_adt)
    '''
    from sklearn.neighbors import NeighborhoodComponentsAnalysis
    adata = adata
    if n_top_genes is not None:
        sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=n_top_genes)
    else:
        if add_features is not None:  # first add the features, input should be anndata
            adata = concat_rna_and_other(adata,add_features,umap=None,name='add_features',prefix='add_features_')
        if hv_features is not None:    # custom hv
            tmp = pd.Series(index=adata.var_names,data=np.full(len(adata.var_names),fill_value=False))
            tmp.loc[hv_features] = True
            adata.var['highly_variable'] = tmp.values
    adata.raw = adata
    adata = adata[:,adata.var['highly_variable']]
    X = make_sure_mat_dense(adata.X)
    y = adata.obs[label].values
    nca = NeighborhoodComponentsAnalysis(n_components=nca_n_components,max_iter=max_iter)
    embed = nca.fit_transform(X,y)  # (n_cells,n_components)
    adata.obsm['X_nca'] = embed
    adata = adata.raw.to_adata()
    if method == 'umap':
        sc.pp.neighbors(adata,use_rep='X_nca')
        sc.tl.umap(adata)
        sc.pl.umap(adata,color=label,frameon=False,legend_loc=legend_loc)
        if save:
            plt.savefig(os.path.join('.','nca_embedding_{}_{}.{}'.format(label,method,format)),bbox_inches='tight')
            plt.close()
    elif method == 'tsne':
        sc.tl.tsne(adata,use_rep='X_nca')
        sc.pl.tsne(adata,color=label,frameon=False,legend_loc=legend_loc)
        if save:
            plt.savefig(os.path.join('.','nca_embedding_{}_{}.{}'.format(label,method,format)),bbox_inches='tight')
            plt.close()
    adata.X = make_sure_mat_sparse(adata.X)
    return adata



def umap_dual_view_save(adata,cols,method='umap'):
    '''
    generate a pdf file with two umap up and down, one is with legend on side, another is with legend on data.
    More importantly, this allows you to generate multiple columns iteratively.

    :param adata: Anndata
    :param cols: list, all columns from which we want to draw umap.
    :param method: string, either umap or tsne

    Examples::

        from sctriangulate.preprocessing import umap_dual_view_save
        umap_dual_view_save(adata,cols=['annotation1','annotation2','total_counts'])
    '''
    if method == 'umap':
        for col in cols:
            fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
            sc.pl.umap(adata,color=col,frameon=False,ax=ax[0])
            sc.pl.umap(adata,color=col,frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
            plt.savefig('./umap_dual_view_{}.pdf'.format(col),bbox_inches='tight')
            plt.close()
    elif method == 'tsne':
        for col in cols:
            fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
            sc.pl.tsne(adata,color=col,frameon=False,ax=ax[0])
            sc.pl.tsne(adata,color=col,frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
            plt.savefig('./tsne_dual_view_{}.pdf'.format(col),bbox_inches='tight')
            plt.close()      

def just_log_norm(adata):
    '''
    perform CPTT and log operation on adata.X in place, no return value

    :param adata: the anndata to perform operations on

    Example:

        just_log_norm(adata)
    '''
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def format_find_concat(adata,canonical_chr_only=True,gtf_file='gencode.v38.annotation.gtf',key_added='gene_annotation',**kwargs):
    '''
    this is a wrapper function to add nearest genes to your ATAC peaks or bins. For instance, if the peak is chr1:55555-55566,
    it will be annotated as chr1:55555-55566_gene1;gene2
    :param adata: The anndata, the var_names is the peak/bin, please make sure the format is like chr1:55555-55566
    :param canonical_chr_only: boolean, default to True, means only contain features on canonical chromosomes. for human, it is chr1-22 and X,Y
    :param gtf_file: the path to the gtf files, we provide the hg38 on this `google drive link <https://drive.google.com/file/d/11gbJl2-wZr3LbpWaU9RiUAGPebqWYi1z/view?usp=sharing>`_ to download
    :param key_added: string, the column name where the gene annotation will be inserted to adata.var, default is 'gene_annotation'
    :return adata: Anndata, the gene annotation will be added to var, and the var_name will be suffixed with gene annotation, if canonical_chr_only is True, then only features on canonical 
                   chromsome will be retained.
    Example::
        adata = format_find_concat(adata)
    '''
    adata= reformat_peak(adata,canonical_chr_only=canonical_chr_only)
    find_genes(adata,gtf_file=gtf_file,key_added=key_added,**kwargs)
    adata.var_names = [name + '_' + gene for name,gene in zip(adata.var_names,adata.var[key_added])]
    return adata



class GeneConvert(object):
    '''
    A collection of gene symbol conversion functions.

    Now support:

    1. ensemblgene id to gene symbol.


    '''

    @staticmethod
    def ensemblgene_to_symbol(query,species):
        '''
        Examples::

            from sctriangulate.preprocessing import GeneConvert
            converted_list = GeneConvert.ensemblgene_to_symbol(['ENSG00000010404','ENSG00000010505'],species='human')
        '''
        # assume query is a list, will also return a list
        import mygene
        mg = mygene.MyGeneInfo()
        out = mg.querymany(query,scopes='ensemblgene',fileds='symbol',species=species,returnall=True,as_dataframe=True,df_index=True)
        result = out['out']['symbol'].fillna('unknown_gene').tolist()
        try:
            assert len(query) == len(result)
        except AssertionError:    # have duplicate results
            df = out['out']
            df_unique = df.loc[~df.index.duplicated(),:]
            result = df_unique['symbol'].fillna('unknown_gene').tolist()
        return result


def dual_gene_plot(adata,gene1,gene2,s=8,save=True,format='pdf',dir='.',umap_lim=None):
    from scipy.sparse import issparse
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    index1 = np.where(adata.var_names == gene1)[0][0]
    index2 = np.where(adata.var_names == gene2)[0][0]
    exp1 = adata.X[:,index1]
    exp2 = adata.X[:,index2]
    color = []
    for i in range(len(exp1)):
        if exp1[i] > 0 and exp2[i] > 0:
            color.append('#F2DE77')
        elif exp1[i] > 0 and exp2[i] == 0:
            color.append('#5ABF9A')
        elif exp1[i] == 0 and exp2[i] > 0:
            color.append('#F25C69')
        else:
            color.append('lightgrey')
    fig, ax = plt.subplots()
    if umap_lim is not None:
        ax.set_xlim(umap_lim[0])
        ax.set_ylim(umap_lim[1])
    ax.scatter(x=adata.obsm['X_umap'][:,0],y=adata.obsm['X_umap'][:,1],s=s,c=color)
    import matplotlib.lines as mlines
    ax.legend(handles=[mlines.Line2D([],[],marker='o',color=i,linestyle='') for i in ['#F2DE77','#5ABF9A','#F25C69','lightgrey']],
              labels=['Both','{}'.format(gene1),'{}'.format(gene2),'None'],frameon=False,loc='upper left',bbox_to_anchor=[1,1])
    if save:
        plt.savefig(os.path.join(dir,'sctri_dual_gene_plot_{}_{}.{}'.format(gene1,gene2,format)),bbox_inches='tight')
        plt.close()
    return ax


def multi_gene_plot(adata,genes,s=8,save=True,format='pdf',dir='.',umap_lim=None):
    from scipy.sparse import issparse
    if issparse(adata.X):
        adata.X = adata.X.toarray()
    exp_list = []
    for gene in genes:
        index_gene = np.where(adata.var_names == gene)[0][0]
        exp_gene = adata.X[:,index_gene]
        exp_list.append(exp_gene)
    color = []
    for i in range(len(exp_list[0])):
        if len(genes) == 3:
            c = ['#04BFBF','#83A603','#F7766D']
        elif len(genes) == 4:
            c = ['#04BFBF','#83A603','#F7766D','#E36DF2']
        elif len(genes) == 5:
            c = ['#04BFBF','#83A603','#F7766D','#E36DF2','#A69B03']                       
        b = '#BABABA'
        l_exp = np.array([exp[i] for exp in exp_list])
        n_exp = np.count_nonzero(l_exp > 0)
        if n_exp > 1:
            color.append(c[np.where(l_exp==l_exp.max())[0][0]])
        elif n_exp == 1:
            color.append(c[np.where(l_exp>0)[0][0]])
        elif n_exp == 0:
            color.append(b)
    fig, ax = plt.subplots()
    if umap_lim is not None:
        ax.set_xlim(umap_lim[0])
        ax.set_ylim(umap_lim[1])
    ax.scatter(x=adata.obsm['X_umap'][:,0],y=adata.obsm['X_umap'][:,1],s=s,c=color)
    import matplotlib.lines as mlines
    ax.legend(handles=[mlines.Line2D([],[],marker='o',color=i,linestyle='') for i in c+[b]],
              labels=genes + ['None'],frameon=False,
              loc='upper left',bbox_to_anchor=[1,1])
    if save:
        output = '_'.join(genes)
        plt.savefig(os.path.join(dir,'sctri_multi_gene_plot_{}.{}'.format(output,format)),bbox_inches='tight')
        plt.close()
    return ax


def make_sure_mat_dense(mat):
    '''
    make sure a matrix is dense

    :param mat: ndarary

    :return mat: ndarray (dense)

    Examples::

        mat = make_sure_mat_dense(mat)
    '''
    if not issparse(mat):
        pass
    else:
        mat = mat.toarray()
    return mat

def make_sure_mat_sparse(mat):  # will be csr if the input mat is a dense array
    '''
    make sure a matrix is sparse

    :param mat: ndarary

    :return mat: ndarray (sparse)

    Examples::

        mat = make_sure_mat_dense(mat)
    '''
    if not issparse(mat):
        mat = csr_matrix(mat)
    else:
        pass
    return mat

class Normalization(object):
    '''
    a series of Normalization functions

    Now support:

    1. CLR normalization
    2. total count normalization (CPTT, CPM)
    3. GMM normalization



    '''
    # matrix should be cell x feature, expecting a ndarray

    @staticmethod
    def CLR_normalization(mat):
        '''
        Examples::

            from sctriangulate.preprocessing import Normalization
            post_mat = Normalization.CLR_normalization(pre_mat)
        '''
        from scipy.stats import gmean
        gmeans = gmean(mat+1,axis=1).reshape(-1,1)
        post = np.log(mat/gmeans + 1)
        return post

    @staticmethod
    def total_normalization(mat,target=1e4):
        '''
        Examples::

            from sctriangulate.preprocessing import Normalization
            post_mat = Normalization.total_normalization(pre_mat)
        '''
        if target is None:
            target = np.mean(mat,axis=1).reshape(-1,1)
        total = np.sum(mat,axis=1).reshape(-1,1)
        sf = total/target
        post = np.log(mat/sf + 1)
        return post

    @staticmethod
    def GMM_normalization(mat,non_negative=False):
        '''
        This method is a re-implementaion from `Stephenson et al <https://www.nature.com/articles/s41591-021-01329-2>`_,
        The raw counts are first subjected to a CPTT total normalization, then a GaussianMixture model was fitted to the data,
        we substract the mean of the background from the data to remove background noise. Optionally, user can make the post-processed
        matrix as non-negative by setting ``non_negative==True``

        Examples::

            from sctriangulate.preprocessing import Normalization
            post_mat = Normalization.GMM_normalization(pre_mat)
        '''
        mat = Normalization.total_normalization(mat)
        from sklearn.mixture import GaussianMixture
        model = GaussianMixture(n_components=2,random_state=0)
        model.fit(mat)
        means = model.means_  # (n_components,n_features)
        bg_index = np.argmin(means.mean(axis=1))
        bg_mean = means[bg_index,:].reshape(1,-1)
        post = mat - bg_mean
        if non_negative:
            post = np.where(post>=0,post,0)
        return post


def gene_activity_count_matrix_new_10x(fall_in_promoter,fall_in_gene,valid=None):
    '''
    Full explanation please refer to ``gene_activity_count_matrix_old_10x``

    Examples::

        from sctriangulate.preprocessing import gene_activity_count_matrix_new_10x
        gene_activity_count_matrix_new_10x(fall_in_promoter,fall_in_gene,valid=None)       
    '''

    gene_promoter = pd.read_csv(fall_in_promoter,sep='\t',header=None)
    gene_body = pd.read_csv(fall_in_gene,sep='\t',header=None)
    bucket = []
    for i in range(gene_promoter.shape[0]):
        row = gene_promoter.iloc[i]
        in_gene = row[3]
        in_barcode = row[6]
        in_count = row[7]
        try:
            in_barcode = in_barcode.split(';')
            in_count = [int(item) for item in in_count.split(';')]
        except AttributeError:  # means no fragments fall into the promoter
            continue
        # tmp will be three column, barcode, count, gene, no index
        tmp = pd.DataFrame({'barcode':in_barcode,'count':in_count}).groupby(by='barcode')['count'].sum().to_frame()
        tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
        tmp.reset_index(inplace=True)
        bucket.append(tmp)
    for i in range(gene_body.shape[0]):
        row = gene_body.iloc[i]
        in_gene = row[3]
        in_barcode = row[6]
        in_count = row[7]
        try:
            in_barcode = in_barcode.split(';')
            in_count = [int(item) for item in in_count.split(';')]
        except AttributeError:  # means no fragments fall into the promoter
            continue
        # tmp will be three column, barcode, count, gene, no index
        tmp = pd.DataFrame({'barcode':in_barcode,'count':in_count}).groupby(by='barcode')['count'].sum().to_frame()
        tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
        tmp.reset_index(inplace=True)
        bucket.append(tmp)
    df = pd.concat(bucket)
    if valid is not None:
        df = df.loc[df['barcode'].isin(valid),:]
    final = df.groupby(by=['barcode','gene'])['count'].sum().unstack(fill_value=0)
    return final

    

def gene_activity_count_matrix_old_10x(fall_in_promoter,fall_in_gene,valid=None):
    '''
    this function is to generate gene activity count matrix, please refer to ``gene_activity_count_matrix_new_10x`` for latest
    version of 10x fragements.tsv output.

    how to get these two arguments? (LIGER team approach)


    1. sort the fragment, gene and promoter bed, or use function in this module to sort the reference bed files::

        sort -k1,1 -k2,2n -k3,3n pbmc_granulocyte_sorted_10k_atac_fragments.tsv > atac_fragments.sort.bed
        sort -k 1,1 -k2,2n -k3,3n hg19_genes.bed > hg19_genes.sort.bed
        sort -k 1,1 -k2,2n -k3,3n hg19_promoters.bed > hg19_promoters.sort.bed

    2. bedmap::

        module load bedops
        bedmap --ec --delim "\t" --echo --echo-map-id hg19_promoters.sort.bed atac_fragments.sort.bed > atac_promoters_bc.bed
        bedmap --ec --delim "\t" --echo --echo-map-id hg19_genes.sort.bed atac_fragments.sort.bed > atac_genes_bc.bed


    the following was taken from http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html

    * **delim**. This changes output delimiter from ‘|’ to indicated delimiter between columns, which in our case is “\t”.

    * **ec**. Adding this will check all problematic input files.

    * **echo**. Adding this will print each line from reference file in output. The reference file in our case is gene or promoter index.

    * **echo-map-id**. Adding this will list IDs of all overlapping elements from mapping files, which in our case are cell barcodes from fragment files.

    3. Finally::

        from sctriangulate.preprocessing import gene_activity_count_matrix_old_10x
        gene_activity_count_matrix_old_10x(fall_in_promoter,fall_in_gene,valid=None)
   
    '''
    gene_promoter = pd.read_csv(fall_in_promoter,sep='\t',header=None)
    gene_body = pd.read_csv(fall_in_gene,sep='\t',header=None)
    bucket = []
    for i in range(gene_promoter.shape[0]):
        row = gene_promoter.iloc[i]
        in_gene = row[3]
        in_barcode = row[6]
        try:
            in_barcode = in_barcode.split(';')
        except AttributeError:  # means no fragments fall into the promoter
            continue
        tmp = pd.Series(in_barcode).value_counts().to_frame(name='count')
        tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
        tmp.reset_index(inplace=True)  # three column: index, count, gene
        bucket.append(tmp)
    for i in range(gene_body.shape[0]):
        row = gene_body.iloc[i]
        in_gene = row[3]
        in_barcode = row[6]
        try:
            in_barcode = in_barcode.split(';')
        except AttributeError:  # means no fragments fall into the promoter
            continue
        tmp = pd.Series(in_barcode).value_counts().to_frame(name='count')
        tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
        tmp.reset_index(inplace=True)  # three column: index, count, gene
        bucket.append(tmp)
    df = pd.concat(bucket)
    if valid is not None:
        df = df.loc[df['index'].isin(valid),:]
    final = df.groupby(by=['index','gene'])['count'].sum().unstack(fill_value=0)
    return final


def gene_bed_to_promoter_bed(gene_bed_path,promoter_bed_path,up_bp=3000):
    gene_bed = pd.read_csv(gene_bed_path,header=None,sep='\t')
    with open(promoter_bed_path,'w') as f:
        for i in range(gene_bed.shape[0]):
            row = gene_bed.iloc[i]
            chro = row[0]
            start = row[1]
            end = row[2]
            name = row[3]
            score = row[4]
            strand = row[5]
            if strand == '+':
                new_start = start - up_bp
                new_end = start
            else:
                new_start = end
                new_end = end + up_bp
            f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chro,new_start,new_end,name,score,strand))


def ensembl_gtf_to_gene_bed(gtf_path,bed_path,sort=True):

    # only interested in gene feature
    hg38_gtf = pd.read_csv(gtf_path,skiprows=5,header=None,sep='\t')
    hg38_gtf_gene = hg38_gtf.loc[hg38_gtf[2]=='gene',:]

    # gotta have gene symbol
    col = []
    for i in range(hg38_gtf_gene.shape[0]):
        metadata = hg38_gtf_gene.iloc[i,8]
        if 'gene_name' in metadata:
            col.append(True)
        else:
            col.append(False)
    hg38_gtf_gene_have_symbol = hg38_gtf_gene.loc[col,:]

    # add biotype and gene name
    col1 = []
    col2 = []
    for i in range(hg38_gtf_gene_have_symbol.shape[0]):
        metadata = hg38_gtf_gene_have_symbol.iloc[i, 8]
        biotype = metadata.split('; ')[-1].split(' ')[-1].strip(';').strip('"')
        name = metadata.split('; ')[2].split(' ')[1].strip('"')
        col1.append(biotype)
        col2.append(name)
    hg38_gtf_gene_have_symbol['biotype'] = col1
    hg38_gtf_gene_have_symbol['name'] = col2

    # biotype has to be either protein_coding, IG or TR gene
    col = (hg38_gtf_gene_have_symbol['biotype']=='protein_coding') |\
          (hg38_gtf_gene_have_symbol['biotype']=='IG_C_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='IG_D_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='IG_J_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='IG_V_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='TR_C_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='TR_D_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='TR_J_gene') |\
          (hg38_gtf_gene_have_symbol['biotype']=='TR_V_gene')
    hg38_gtf_gene_have_symbol_biotype_correct = hg38_gtf_gene_have_symbol.loc[col,:]

    # chromsome need to be correct
    chr_need_chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
    chr_need_int = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    chr_need_other = ['X','Y'] # don't include MT decause fragment.tsv file doesn't output that
    chr_need = chr_need_chr + chr_need_int + chr_need_other

    hg38_gtf_gene_have_symbol_biotype_correct_chr = hg38_gtf_gene_have_symbol_biotype_correct.loc[hg38_gtf_gene_have_symbol_biotype_correct[0].isin(chr_need),:]

    prefixed_chr = ['chr' + str(item) for item in hg38_gtf_gene_have_symbol_biotype_correct_chr[0]]
    hg38_gtf_gene_have_symbol_biotype_correct_chr[0] = prefixed_chr

    # get final result, BED6 format
    final = hg38_gtf_gene_have_symbol_biotype_correct_chr.loc[:,[0,3,4,'name',5,6]]
    if sort:
        '''
        equivalent to:
        sort -k1,1 -k2,2n -k3,3n gene.bed
        '''
        final.sort_values(by=[0,3,4],inplace=True)
    final.to_csv(bed_path,sep='\t',index=None,header=None)

# this function is taken from episcanpy, all the credits to the original developer:
# https://github.com/colomemaria/epiScanpy/blob/master/episcanpy/tools/_find_genes.py

# to download the gtf file
'''
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
'''
def find_genes(adata,
                 gtf_file,
                 key_added='gene_annotation',
                 upstream=2000,
                 downstream=0,
                 feature_type='gene',
                 annotation='HAVANA',
                 raw=False):
    """
    This function is taken from `episcanpy <https://github.com/colomemaria/epiScanpy/blob/master/episcanpy/tools/_find_genes.py>`_.
    all the credits to the original developer

    merge values of peaks/windows/features overlapping genebodies + 2kb upstream. It is possible to extend the search for closest gene to a given 
    number of bases downstream as well. There is commonly 2 set of annotations in a gtf file(HAVANA, ENSEMBL). By default, the function will search 
    annotation from HAVANA but other annotation label/source can be specifed. It is possible to use other type of features than genes present in a 
    gtf file such as transcripts or CDS.


    Examples::

        from sctriangulate.preprocessing import find_genes
        find_genes(adata,gtf_file='gencode.v38.annotation.gtf)
    
    """
    ### extracting the genes
    gtf = {}
    with open(gtf_file) as f:
        for line in f:
            if line[0:2] != '##' and '\t'+feature_type+'\t' in line and '\t'+annotation+'\t' in line:
                line = line.rstrip('\n').split('\t')
                if line[6] == '-':
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3])-downstream, int(line[4])+upstream,line[-1].split(';')[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3])-downstream, int(line[4])+upstream,line[-1].split(';')[:-1]])
                else:
                    if line[0] not in gtf.keys():
                        gtf[line[0]] = [[int(line[3])-upstream, int(line[4])+downstream,line[-1].split(';')[:-1]]]
                    else:
                        gtf[line[0]].append([int(line[3])-upstream, int(line[4])+downstream,line[-1].split(';')[:-1]])

    # extracting the feature coordinates
    raw_adata_features = {}
    feature_index = 0
    for line in adata.var_names.tolist():
        line = line.split('_')
        if line[0] not in raw_adata_features.keys():
            raw_adata_features[line[0]] = [[int(line[1]),int(line[2]), feature_index]]
        else:
            raw_adata_features[line[0]].append([int(line[1]),int(line[2]), feature_index])
        feature_index += 1
    
    ## find the genes overlaping the features.
    gene_index = []
    for chrom in raw_adata_features.keys():
        if chrom in gtf.keys():
            chrom_index = 0
            previous_features_index = 0
            for feature in raw_adata_features[chrom]:
                gene_name = []
                feature_start = feature[0]
                feature_end = feature[1]
                for gene in gtf[chrom]:
                    if (gene[1]<= feature_start): # the gene is before the feature. we need to test the next gene.
                        continue
                    elif (feature_end <= gene[0]): # the gene is after the feature. we need to test the next feature.
                        break
                    else: # the window is overlapping the gene. 
                        for n in gene[-1]:
                            if 'gene_name' in n:
                                gene_name.append(n.lstrip('gene_name "').rstrip('""'))
                        
                if gene_name == []:
                    gene_index.append('intergenic')
                elif len(gene_name)==1:
                    gene_index.append(gene_name[0])
                else:
                    gene_index.append(";".join(list(set(gene_name))))
                    
    adata.var[key_added] = gene_index

    
def reformat_peak(adata,canonical_chr_only=True):
    '''
    To use ``find_genes`` function, please first reformat the peak from 10X format "chr1:10109-10357" to
    find_gene format "chr1_10109_10357"

    :param adata: AnnData
    :param canonical_chr_only: boolean, only kept the canonical chromosome

    :return: AnnData

    Examples::

        from sctriangulate.preprocessing import reformat_peak
        adata = reformat_peak(adata,canonical_chr_only=True)

    '''
    var_names = adata.var_names
    valid = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16',
             'chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])
    col = []
    for item in var_names:
        chr_ = item.split(':')[0]
        if chr_ in valid:
            start = item.split(':')[1].split('-')[0]
            end = item.split(':')[1].split('-')[1]
            now = '_'.join([chr_,start,end])
            col.append(now)
        else:
            col.append('non_canonical_chr')
    adata.var_names = col
    if canonical_chr_only:
        adata = adata[:,adata.var_names!='non_canonical_chr'].copy()
    return adata


def plot_coexpression(adata,gene1,gene2,kind,hist2d_bins=50,hist2d_cmap=bg_greyed_cmap('viridis'),hist2d_vmin=1e-5,hist2d_vmax=None,
                      scatter_dot_color='blue',contour_cmap='viridis',contour_levels=None,contour_scatter=True,contour_scatter_dot_size=5,
                      contour_train_kde='valid',surface3d_cmap='coolwarm',save=True,outdir='.',name=None):
    x = np.squeeze(make_sure_mat_dense(adata[:,gene1].X))
    y = np.squeeze(make_sure_mat_dense(adata[:,gene2].X))
    if kind == 'scatter':
        fig,ax = plt.subplots()
        ax.scatter(x,y,color=scatter_dot_color)
        ax.set_xlabel('{}'.format(gene1))
        ax.set_ylabel('{}'.format(gene2))
    elif kind == 'hist2d':
        fig,ax = plt.subplots()
        hist2d = ax.hist2d(x,y,bins=hist2d_bins,cmap=hist2d_cmap,vmin=hist2d_vmin,vmax=hist2d_vmax)
        fig.colorbar(mappable=hist2d[3],ax=ax)
        ax.set_xlabel('{}'.format(gene1))
        ax.set_ylabel('{}'.format(gene2))
    elif kind == 'contour':
        from scipy.stats import gaussian_kde
        fig,ax = plt.subplots()
        X,Y = np.meshgrid(np.linspace(x.min(),x.max(),100),np.linspace(y.min(),y.max(),100))
        positions = np.vstack([X.ravel(),Y.ravel()])  # (2, 10000)
        values = np.vstack([x,y])   # (2, 2700)
        if contour_train_kde == 'valid':  # data points that are non-zero for both gene1 and gen2
            values_to_kde = values[:,np.logical_not(np.any(values==0,axis=0))]
        elif contour_train_kde == 'semi_valid':  # data points that are non-zero for at least one of the gene
            values_to_kde = values[:,np.logical_not(np.all(values==0,axis=0))]
        elif contour_train_kde == 'full':   # all data points will be used for kde estimation
            values_to_kde == values
        kernel = gaussian_kde(values_to_kde)  
        density = kernel(positions) # (10000,)
        density = density.reshape(X.shape)  # (100,100)
        cset = ax.contour(X,Y,density,levels=contour_levels,cmap=contour_cmap)
        if contour_scatter:
            dot_density = kernel(values)
            dot_density_color = [cm.viridis(round(np.interp(x=item,xp=[dot_density.min(),dot_density.max()],fp=[0,255]))) for item in dot_density]
            ax.scatter(x,y,c=dot_density_color,s=contour_scatter_dot_size)
        from matplotlib.colors import Normalize
        fig.colorbar(mappable=cm.ScalarMappable(norm=Normalize(),cmap=contour_cmap),ax=ax)
        ax.set_xlabel('{}'.format(gene1))
        ax.set_ylabel('{}'.format(gene2))

    elif kind == 'contourf':
        from scipy.stats import gaussian_kde
        fig,ax = plt.subplots()
        X,Y = np.meshgrid(np.linspace(x.min(),x.max(),100),np.linspace(y.min(),y.max(),100))
        positions = np.vstack([X.ravel(),Y.ravel()])  # (2, 10000)
        values = np.vstack([x,y])   # (2, 2700)
        if contour_train_kde == 'valid':  # data points that are non-zero for both gene1 and gen2
            values_to_kde = values[:,np.logical_not(np.any(values==0,axis=0))]
        elif contour_train_kde == 'semi_valid':  # data points that are non-zero for at least one of the gene
            values_to_kde = values[:,np.logical_not(np.all(values==0,axis=0))]
        elif contour_train_kde == 'full':   # all data points will be used for kde estimation
            values_to_kde == values
        kernel = gaussian_kde(values_to_kde)  
        density = kernel(positions) # (10000,)
        density = density.reshape(X.shape)  # (100,100)
        cfset = ax.contourf(X,Y,density,levels=contour_levels,cmap=contour_cmap) 
        cset = ax.contour(X,Y,density,levels=contour_levels,colors='k')  
        clable = ax.clabel(cset,inline=True,fontsize=5)
        from matplotlib.colors import Normalize
        fig.colorbar(mappable=cm.ScalarMappable(norm=Normalize(),cmap=contour_cmap),ax=ax)
        ax.set_xlabel('{}'.format(gene1))
        ax.set_ylabel('{}'.format(gene2))

    elif kind == 'surface3d':
        fig = plt.figure()
        from scipy.stats import gaussian_kde
        X,Y = np.meshgrid(np.linspace(x.min(),x.max(),100),np.linspace(y.min(),y.max(),100))
        positions = np.vstack([X.ravel(),Y.ravel()])  # (2, 10000)
        values = np.vstack([x,y])   # (2, 2700)
        if contour_train_kde == 'valid':  # data points that are non-zero for both gene1 and gen2
            values_to_kde = values[:,np.logical_not(np.any(values==0,axis=0))]
        elif contour_train_kde == 'semi_valid':  # data points that are non-zero for at least one of the gene
            values_to_kde = values[:,np.logical_not(np.all(values==0,axis=0))]
        elif contour_train_kde == 'full':   # all data points will be used for kde estimation
            values_to_kde == values
        kernel = gaussian_kde(values_to_kde)  
        density = kernel(positions) # (10000,)
        density = density.reshape(X.shape)  # (100,100)
        ax = plt.axes(projection='3d')
        surf = ax.plot_surface(X,Y,density,cmap=surface3d_cmap)
        ax.set_xlabel('{}'.format(gene1))
        ax.set_ylabel('{}'.format(gene2))  
        ax.set_zlabel('PDF for KDE')      
        fig.colorbar(mappable=surf,ax=ax)
    if save:
        if name is None:
            plt.savefig(os.path.join(outdir,'coexpression_{}_{}_{}_plot.pdf'.format(kind,gene1,gene2)),bbox_inches='tight')
            plt.close()
        else:
            plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
            plt.close()            

    return ax


def umap_color_exceed_102(adata,key,dot_size=None,legend_fontsize=6,outdir='.',name=None):
    '''
    draw a umap that bypass the scanpy 102 color upper bound, this can generate as many as 433 clusters.
    :param adata: Anndata
    :param key: the categorical column in adata.obs, which will be plotted
    :param dot_size: None or number
    :param legend_fontsize: defualt is 6
    :param outdir: output directory, default is '.'
    :param name: name of the plot, default is None

    Exmaple::

        from sctriangulate.preprocessing import umap_color_exceed_102
        umap_color_exceed_102(adata,key='leiden6')  # more than 130 clusters

    .. image:: ./_static/more_than102.png
        :height: 550px
        :width: 550px
        :align: center
        :target: target       
    '''
    fig,ax = plt.subplots()
    mapping = colors_for_set(adata.obs[key].unique().tolist())
    color = adata.obs[key].map(mapping).values
    if dot_size is None:
        dot_size = 120000/adata.shape[0]
    ax.scatter(adata.obsm['X_umap'][:,0],adata.obsm['X_umap'][:,1],c=color,s=dot_size)
    import matplotlib.lines as mlines
    ax.legend(handles=[mlines.Line2D([],[],marker='o',linestyle='',color=i) for i in mapping.values()],
            labels=[i for i in mapping.keys()],loc='upper left',bbox_to_anchor=(1,1),ncol=3,frameon=False,prop={'size':6})
    if name is None:
        name = 'umap_{}_exceed_102.pdf'.format(key)
    plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
    plt.close()


def custom_two_column_sankey(adata,left_annotation,right_annotation,opacity=0.6,pad=3,thickness=10,margin=300,text=True,save=True,as_html=True,outdir='.'):
    import plotly.graph_objects as go
    import kaleido
    df = adata.obs.loc[:,[left_annotation,right_annotation]]
    node_label = df[left_annotation].unique().tolist() + df[right_annotation].unique().tolist()
    node_color = pick_n_colors(len(node_label))
    link = []
    for source,sub in df.groupby(by=left_annotation):
        for target,subsub in sub.groupby(by=right_annotation):
            if subsub.shape[0] > 0:
                link.append((source,target,subsub.shape[0],subsub.shape[0]/sub.shape[0]))
    link_info = list(zip(*link))
    link_source = [node_label.index(item) for item in link_info[0]]
    link_target = [node_label.index(item) for item in link_info[1]]
    link_value = link_info[2]
    link_pert = [round(item,2) for item in link_info[3]]
    
    link_color = ['rgba{}'.format(tuple([infer_to_256(item) for item in to_rgb(node_color[i])] + [opacity])) for i in link_source]
    node_plotly = dict(pad = pad, thickness = thickness,line = dict(color = "grey", width = 0.1),label = node_label,color = node_color)
    link_plotly = dict(source=link_source,target=link_target,value=link_value,color=link_color,customdata=link_pert,
                       hovertemplate='%{source.label} -> %{target.label} <br /> number of cells: %{value} <br /> percentage: %{customdata}')
    if not text:
        fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly, textfont=dict(color='rgba(0,0,0,0)',size=1))])
    else:
        fig = go.Figure(data=[go.Sankey(node = node_plotly,link = link_plotly)])
    fig.update_layout(title_text='sankey_{}_{}'.format(left_annotation,right_annotation), font_size=6, margin=dict(l=margin,r=margin))
    if save:
        if not as_html:
            fig.write_image(os.path.join(outdir,'two_column_sankey_{}_{}_text_{}.pdf'.format(left_annotation,right_annotation,text))) 
        else:
            fig.write_html(os.path.join(outdir,'two_column_sankey_{}_{}_text_{}.html'.format(left_annotation,right_annotation,text)),include_plotlyjs='cdn') 


def rna_umap_transform(outdir,ref_exp,ref_group,q_exp_list,q_group_list,q_identifier_list,pca_n_components,umap_min_dist=0.5):
    '''
    Take a reference expression matrix (pandas dataframe), and a list of query expression matrics (pandas dataframe),
    along with a list of query expression identifiers. This function will generate the umap transformation of your reference exp,
    and make sure to squeeze your each query exp to the same umap transformation.

    :param outdir: the path in which all the results will go
    :param ref_exp: the pandas dataframe object,index is the features, column is the cell barcodes
    :param ref_group: the pandas series object, index is the cell barcodes, the value is the clusterlabel information
    :param q_exp_list: the list of pandas dataframe object, requirement is the same as ref_exp
    :param q_group_list: the list of pandas series object, requirement is the same as ref_group
    :param q_identifier_list: the list of string, to denote the name of each query dataset
    :param pca_n_components: the int, denote the PCs to use for PCA
    :param umap_min_dist: the float number, denote the min_dist parameter for umap program

    Examples::

        rna_umap_transform(outdir='.',ref_exp=ref_df,ref_group=ref_group_df,q_exp_list=[q1_df],q_group_list=[q1_group_df],q_identifier_list=['q1'],pca_n_components=50)
    '''
    # create outdir if not exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ref_variable = ref_exp.index.values
    store_df = []
    store_variable = []
    for q_exp in q_exp_list:
        q_variable = q_exp.index.values
        store_df.append(q_exp)
        store_variable.append(q_variable)

    # find common features between ref and all qs.
    all_variable = copy.deepcopy(store_variable)
    all_variable.insert(0,ref_variable)
    common_variable = list(reduce(lambda a,b: set(a).intersection(set(b)),all_variable))

    # subset the exps, also transpose
    ref_exp_common = ref_exp.loc[common_variable,:].T
    store_df_common = []
    for q_exp in store_df:
        q_exp_common = q_exp.loc[common_variable,:].T
        store_df_common.append(q_exp_common)


    # train pca and umap on ref
    ref_pca_model = PCA(n_components = pca_n_components).fit(ref_exp_common.values)
    ref_pca_score = ref_pca_model.transform(ref_exp_common.values)
    ref_umap_model = umap.UMAP(min_dist=umap_min_dist).fit(ref_pca_score)
    ref_umap_embed = ref_umap_model.embedding_

    # transform all the querys
    store_q_umap = []
    for q_exp_common in store_df_common:
        q_pca_score = ref_pca_model.transform(q_exp_common.values)
        q_umap_embed = ref_umap_model.transform(q_pca_score)
        store_q_umap.append(q_umap_embed)


    ref_umap_embed_df = pd.DataFrame(data=ref_umap_embed,index=ref_exp_common.index,columns=['umap_x','umap_y'])
    ref_umap_embed_df.to_csv(os.path.join(outdir,'ref_umap.txt'),sep='\t')
    for i,item in enumerate(store_q_umap):
        q_exp_common = store_df_common[i]
        item = pd.DataFrame(data=item,index=q_exp_common.index,columns=['umap_x','umap_y'])
        item.to_csv(os.path.join(outdir,'query_{}_umap.txt'.format(q_identifier_list[i])),sep='\t')


    # visualization
    all_identifier = ['ref'] + q_identifier_list
    all_exp = [ref_exp_common] + store_df_common
    all_label_mapping = [group_df.to_dict() for group_df in [ref_group] + q_group_list]
    all_umap = [ref_umap_embed] + store_q_umap
    for i,exp,label_map,embed in zip(all_identifier,all_exp,all_label_mapping,all_umap):
        adata = ad.AnnData(X=exp.values,obs=pd.DataFrame(index=exp.index),var=pd.DataFrame(index=exp.columns))
        adata.obs['label'] = adata.obs_names.map(label_map).fillna('unknown').values
        adata.obsm['X_umap'] = embed
        sc.pl.umap(adata,color='label',legend_loc='on data',legend_fontsize=6)
        plt.savefig(os.path.join(outdir,'umap_{}.pdf'.format(i)),bbox_inches='tight')
        plt.close()