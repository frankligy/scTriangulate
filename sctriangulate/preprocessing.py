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


def large_txt_to_mtx(int_file,out_folder,gene_is_index=True):
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


def mtx_to_large_txt(int_folder,out_file,gene_is_index=False):
    gene = pd.read_csv(os.path.join(int_folder,'genes.tsv'),sep='\t',index_col=0,header=None).index
    cell = pd.read_csv(os.path.join(int_folder,'barcodes.tsv'),sep='\t',index_col=0,header=None).index
    value = mmread(os.path.join(int_folder,'matrix.mtx')).toarray()
    if gene_is_index:
        data = pd.DataFrame(data=value,index=gene,columns=cell)
    else:
        data = pd.DataFrame(data=value.T,index=cell,columns=gene)
    data.to_csv(out_file,sep='\t',chunksize=1000)


def add_annotations(adata,inputs,cols=None):
    # means a single file such that first column is barcodes, annotations are in the following columns
    annotations = pd.read_csv(inputs,sep='\t',index_col=0).loc[cols]
    mappings = []
    for col in cols:
        mappping = annotations[col].to_dict()
        mappings.append(mapping)
    for i,col in enumerate(cols):
        adata.obs[col] = adata.obs_names.map(mappings[i]).values
    return adata

def add_umap(adata,inputs,mode,cols=None):
    # make sure cols are [umap_x, umap_y]
    if mode == 'pandas':
        df = pd.read_csv(inputs,sep='\t')
        umap = df[cols].values
        adata.obsm['X_umap'] = umap
    elif mode == 'numpy':
        adata.obsm['X_umap'] = inputs
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





def rna_10x_recipe(h5file,write,output=None,azimuth=None):
    # reading
    adata = sc.read_10x_h5(h5file)
    adata.var_names_make_unique()
    # normal analysis
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
    sc.pp.neighbors(adata,n_pcs=40,n_neighbors=15)
    sc.tl.leiden(adata,resolution=1,key_added='leiden1')
    sc.tl.leiden(adata,resolution=2,key_added='leiden2')
    sc.tl.leiden(adata,resolution=3,key_added='leiden3')
    sc.tl.umap(adata)
    # put azimuth result to adata, using azimuth default options, l2 annotation
    if azimuth is not None:
        azimuth = pd.read_csv(azimuth,sep='\t',index_col=0)
        azimuth_map = azimuth['predicted.celltype.l2'].to_dict()
        azimuth_prediction = azimuth['predicted.celltype.l2.score'].to_dict()
        azimuth_mapping = azimuth['mapping.score'].to_dict()
        adata.obs['azimuth'] = adata.obs_names.map(azimuth_map).values
        adata.obs['prediction_score'] = adata.obs_names.map(azimuth_prediction).values
        adata.obs['mapping_score'] = adata.obs_names.map(azimuth_mapping).values
    adata = adata.raw.to_adata()
    if write:
        adata.write(output)
    return adata

class Normalization(object):
    '''matrix should be cell x ADT, expecting a ndarray'''

    @staticmethod
    def CLR_normalization(mat):
        from scipy.stats import gmean
        gmeans = gmean(mat+1,axis=1).reshape(-1,1)
        post = np.log(mat/gmeans + 1)
        return post

    @staticmethod
    def total_normalization(mat,target=1e6):
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

    




