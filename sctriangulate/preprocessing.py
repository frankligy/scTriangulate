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


def small_txt_to_adata(int_file,gene_is_index=True):
    df = pd.read_csv(int_file,sep='\t',index_col=0)
    if gene_is_index:
        adata = ad.AnnData(X=csr_matrix(df.values.T),var=pd.DataFrame(index=df.index),obs=pd.DataFrame(index=df.columns))
    else:
        adata = ad.AnnData(X=csr_matrix(df.values),var=pd.DataFrame(index=df.columns),obs=pd.DataFrame(index=df.index))
    adata.var_names_make_unique()
    return adata


def large_txt_to_mtx(int_file,out_folder,gene_is_index=True):  # whether the txt if gene * cell
    '''since expression matrix is too large, I need to do iteration'''
    reader = pd.read_csv(int_file,sep='\t',index_col=0,chunksize=1000)
    store = []
    for chunk in reader:
        tmp = chunk.astype('int16')
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


def mtx_to_adata(int_folder,gene_is_index=True):  # whether the mtx file is gene * cell
    gene = pd.read_csv(os.path.join(int_folder,'genes.tsv'),sep='\t',index_col=0,header=None).index
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
    gene = pd.read_csv(os.path.join(int_folder,'genes.tsv'),sep='\t',index_col=0,header=None).index
    cell = pd.read_csv(os.path.join(int_folder,'barcodes.tsv'),sep='\t',index_col=0,header=None).index
    value = mmread(os.path.join(int_folder,'matrix.mtx')).toarray()
    if gene_is_index:
        data = pd.DataFrame(data=value,index=gene,columns=cell)
    else:
        data = pd.DataFrame(data=value.T,index=cell,columns=gene)
    data.to_csv(out_file,sep='\t',chunksize=1000)

def adding_azimuth(adata,result,name='predicted.celltype.l2'):
    azimuth = pd.read_csv(result,sep='\t',index_col=0)
    azimuth_map = azimuth[name].to_dict()
    azimuth_prediction = azimuth['{}.score'.format(name)].to_dict()
    azimuth_mapping = azimuth['mapping.score'].to_dict()
    adata.obs['azimuth'] = adata.obs_names.map(azimuth_map).values
    adata.obs['prediction_score'] = adata.obs_names.map(azimuth_prediction).values
    adata.obs['mapping_score'] = adata.obs_names.map(azimuth_mapping).values

def add_annotations(adata,inputs,cols=None,index_col=0):
    # means a single file such that one column is barcodes, annotations are within other columns
    annotations = pd.read_csv(inputs,sep='\t',index_col=index_col).loc[:,cols]
    mappings = []
    for col in cols:
        mapping = annotations[col].to_dict()
        mappings.append(mapping)
    for i,col in enumerate(cols):
        adata.obs[col] = adata.obs_names.map(mappings[i]).fillna('Unknown').values


def add_umap(adata,inputs,mode,cols=None,index_col=0):
    # make sure cols are [umap_x, umap_y]
    if mode == 'pandas':
        df = pd.read_csv(inputs,sep='\t',index_col=index_col)
        umap_x = df[cols[0]].to_dict()
        umap_y = df[cols[1]].to_dict()
        adata.obs['umap_x'] = adata.obs_names.map(umap_x).values
        adata.obs['umap_y'] = adata.obs_names.map(umap_y).values
        adata.obsm['X_umap'] = adata.obs.loc[:,['umap_x','umap_y']].values
        adata.obs.drop(columns=['umap_x','umap_y'],inplace=True)
    elif mode == 'numpy':  # assume the order is correct
        adata.obsm['X_umap'] = inputs


def scanpy_recipe(adata,is_log,resolutions=[0.5,1,2],modality='rna',umap=True,save=True):
    adata.var_names_make_unique()
    # normal analysis
    if modality == 'rna':
        if not is_log:   # count data
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
            sc.pp.normalize_total(adata,target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
            sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata)
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
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
            sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata)
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
            pass
        else:
            sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
            adata.raw = adata
            adata = adata[:,adata.var['highly_variable']]
            sc.pp.scale(adata,max_value=10)
            sc.tl.pca(adata)
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
            pass
        else:
            sc.tl.pca(adata)
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
    




def concat_rna_and_other(adata_rna,adata_other,umap,name):
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
    adata_combine = ad.concat([adata_rna,adata_other],axis=1,join='outer',merge='first',label='modality',keys=['rna','{}'.format(name)])
    if umap == 'rna':
        adata_combine.obsm['X_umap'] = adata_rna.obsm['X_umap']
    elif umap == 'other':
        adata_combine.obsm['X_umap'] = adata_other.obsm['X_umap']
    return adata_combine






def umap_dual_view_save(adata,cols):
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

    @staticmethod
    def ensemblgene_to_symbol(query,species):
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






def make_sure_mat_dense(mat):
    if not issparse(mat):
        pass
    else:
        mat = mat.toarray()
    return mat

def make_sure_mat_sparse(mat):  # will be csr if the input mat is a dense array
    if not issparse(mat):
        mat = csr_matrix(mat)
    else:
        pass
    return mat

class Normalization(object):
    '''matrix should be cell x feature, expecting a ndarray'''

    @staticmethod
    def CLR_normalization(mat):
        from scipy.stats import gmean
        gmeans = gmean(mat+1,axis=1).reshape(-1,1)
        post = np.log(mat/gmeans + 1)
        return post

    @staticmethod
    def total_normalization(mat,target=1e4):
        total = np.sum(mat,axis=1).reshape(-1,1)
        sf = total/target
        post = np.log(mat/sf + 1)
        return post

    @staticmethod
    def GMM_normalization(mat):
        mat = Normalization.total_normalization(mat)
        from sklearn.mixture import GaussianMixture
        model = GaussianMixture(n_components=2,random_state=0)
        model.fit(mat)
        means = model.means_  # (n_components,n_features)
        bg_index = np.argmin(means.mean(axis=1))
        bg_mean = means[bg_index,:].reshape(1,-1)
        post = mat - bg_mean
        return post
    

def gene_activity_count_matrix(fall_in_promoter,fall_in_gene,valid=None):
    '''
    how to get these two arguments? (LIGER team approach)
    - sort the fragment, gene and promoter bed
    sort -k1,1 -k2,2n -k3,3n pbmc_granulocyte_sorted_10k_atac_fragments.tsv > atac_fragments.sort.bed
    sort -k 1,1 -k2,2n -k3,3n hg19_genes.bed > hg19_genes.sort.bed
    sort -k 1,1 -k2,2n -k3,3n hg19_promoters.bed > hg19_promoters.sort.bed

    - bedmap
    module load bedops
    bedmap --ec --delim "\t" --echo --echo-map-id hg19_promoters.sort.bed atac_fragments.sort.bed > atac_promoters_bc.bed
    bedmap --ec --delim "\t" --echo --echo-map-id hg19_genes.sort.bed atac_fragments.sort.bed > atac_genes_bc.bed
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

    




