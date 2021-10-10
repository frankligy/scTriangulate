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


def mtx_to_adata(int_folder,gene_is_index=True,feature='genes'):  # whether the mtx file is gene * cell
    '''
    convert mtx file to adata in RAM, make sure the X is sparse.

    :param int_folder: string, folder where the mtx files are stored.
    :param gene_is_index: boolean, whether the gene is index.
    :param features: string, the name of the feature tsv file, if rna, it will be genes.tsv.

    :return: AnnData

    Examples::

        from sctriangulate.preprocessing import mtx_to_adata
        mtx_to_adata(int_folder='./data',gene_is_index=False,feature='genes')

    '''
    gene = pd.read_csv(os.path.join(int_folder,'{}.tsv'.format(feature)),sep='\t',index_col=0,header=None).index
    cell = pd.read_csv(os.path.join(int_folder,'barcodes.tsv'),sep='\t',index_col=0,header=None).index
    value = mmread(os.path.join(int_folder,'matrix.mtx')).toarray()
    if gene_is_index:
        data = pd.DataFrame(data=value.T,index=cell,columns=gene)
        adata = ad.AnnData(X=data.values,obs=pd.DataFrame(index=data.index.values),var=pd.DataFrame(index=data.columns.values))
    else:
        data = pd.DataFrame(data=value,index=cell,columns=gene)
        adata = ad.AnnData(X=data.values,obs=pd.DataFrame(index=data.index.values),var=pd.DataFrame(index=data.columns.values))
    adata.X = csr_matrix(adata.X)
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

def add_annotations(adata,inputs,cols_input,index_col=0,cols_output=None):
    '''
    Adding annotations from external sources to the adata

    :param adata: Anndata
    :param inputs: string, path to the txt file where the barcode to cluster label information is stored.
    :param cols_input: list, what columns the users want to transfer to the adata.
    :param index_col: int, for the input, which column will serve as the index column
    :param cols_output: list, corresponding to the cols_input, how these columns will be named in the adata.obs columns

    Examples::

        from sctriangulate.preprocessing import add_annotations
        add_annotations(adata,input='./annotation.txt',cols_input=['col1','col2'],index_col=0,cols_output=['annotation1','annontation2'])

    '''
    # means a single file such that one column is barcodes, annotations are within other columns
    annotations = pd.read_csv(inputs,sep='\t',index_col=index_col).loc[:,cols_input]
    mappings = []
    for col in cols_input:
        mapping = annotations[col].to_dict()
        mappings.append(mapping)
    if cols_output is None:
        for i,col in enumerate(cols_input):
            adata.obs[col] = adata.obs_names.map(mappings[i]).fillna('Unknown').values
    else:
        for i in range(len(cols_input)):
            adata.obs[cols_output[i]] = adata.obs_names.map(mappings[i]).fillna('Unknown').values


def add_umap(adata,inputs,mode,cols=None,index_col=0):
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

    Examples::

        from sctriangulate.preprocessing import add_umap
        add_umap(adata,inputs='umap.txt',mode='pandas_disk',cols=['umap1','umap2],index_col=0)

    '''
    # make sure cols are [umap_x, umap_y]
    if mode == 'pandas_disk':
        df = pd.read_csv(inputs,sep='\t',index_col=index_col)
        umap_x = df[cols[0]].to_dict()
        umap_y = df[cols[1]].to_dict()
        adata.obs['umap_x'] = adata.obs_names.map(umap_x).values
        adata.obs['umap_y'] = adata.obs_names.map(umap_y).values
        adata.obsm['X_umap'] = adata.obs.loc[:,['umap_x','umap_y']].values
        adata.obs.drop(columns=['umap_x','umap_y'],inplace=True)
    elif mode == 'pandas_memory':
        df = inputs
        umap_x = df[cols[0]].to_dict()
        umap_y = df[cols[1]].to_dict()
        adata.obs['umap_x'] = adata.obs_names.map(umap_x).values
        adata.obs['umap_y'] = adata.obs_names.map(umap_y).values
        adata.obsm['X_umap'] = adata.obs.loc[:,['umap_x','umap_y']].values
        adata.obs.drop(columns=['umap_x','umap_y'],inplace=True)
    elif mode == 'numpy':  # assume the order is correct
        adata.obsm['X_umap'] = inputs

def doublet_predict(adata):  # gave RNA count or log matrix
    '''
    wrapper function for running srublet, a new column named 'doublet_scores' will be added to the adata

    :param adata: Anndata

    :return: AnnData

    Examples::

        from sctriangulate.preprocessing import doublet_predict
        new_adata = doublet_predict(old_adata)

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
    maks sure the adata is able to write to disk, since h5 file is stricted typed, so on mixed dtype is allowd.
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


def scanpy_recipe(adata,is_log=False,resolutions=[1,2,3],modality='rna',umap=True,save=True,pca_n_comps=None,n_top_genes=3000):
    '''
    Main preprocessing function. Run Scanpy normal pipeline to achieve Leiden clustering with various resolutions across multiple modalities.

    :param adata: Anndata
    :param is_log: boolean, whether the adata.X is count or normalized data.
    :param resolutions: list, what leiden resolutions the users want to obtain.
    :param modality: string, valid values: 'rna','adt','atac'
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

    '''
    adata.var_names_make_unique()
    # normal analysis
    if modality == 'rna':
        if not is_log:   # count data
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
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
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
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

    return adata
    




def concat_rna_and_other(adata_rna,adata_other,umap,name,prefix):
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
        adata_combine.obsm['X_umap'] = adata_rna.obsm['X_umap']
    elif umap == 'other':
        adata_combine.obsm['X_umap'] = adata_other.obsm['X_umap']
    if not issparse(adata_combine.X):
        adata_combine.X = csr_matrix(adata_combine.X)
    return adata_combine






def umap_dual_view_save(adata,cols):
    '''
    generate a pdf file with two umap up and down, one is with legend on side, another is with legend on data.
    More importantly, this allows you to generate multiple columns iteratively.

    :param adata: Anndata
    :param cols: list, all columns from which we want to draw umap.

    Examples::

        from sctriangulate.preprocessing import umap_dual_view_save
        umap_dual_view_save(adata,cols=['annotation1','annotation2','total_counts'])
    '''
    for col in cols:
        fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
        sc.pl.umap(adata,color=col,frameon=False,ax=ax[0])
        sc.pl.umap(adata,color=col,frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
        plt.savefig('./umap_dual_view_{}.pdf'.format(col),bbox_inches='tight')
        plt.close()


def just_log_norm(adata):
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
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
        total = np.sum(mat,axis=1).reshape(-1,1)
        sf = total/target
        post = np.log(mat/sf + 1)
        return post

    @staticmethod
    def GMM_normalization(mat):
        '''
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






