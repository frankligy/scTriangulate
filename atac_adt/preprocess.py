#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
from pandas.core.groupby.generic import AggScalar
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
from itertools import repeat

# adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)
# valid = adata.obs_names

# gene_promoter = pd.read_csv('atac_promoters_bc.bed',sep='\t',header=None)
# gene_body = pd.read_csv('atac_genes_bc.bed',sep='\t',header=None)
# bucket = []
# for i in range(gene_promoter.shape[0]):
#     row = gene_promoter.iloc[i]
#     in_gene = row[3]
#     in_barcode = row[6]
#     try:
#         in_barcode = in_barcode.split(';')
#     except AttributeError:  # means no fragments fall into the promoter
#         continue
#     tmp = pd.Series(in_barcode).value_counts().to_frame(name='count')
#     tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
#     tmp.reset_index(inplace=True)  # three column: index, count, gene
#     bucket.append(tmp)
# for i in range(gene_body.shape[0]):
#     row = gene_body.iloc[i]
#     in_gene = row[3]
#     in_barcode = row[6]
#     try:
#         in_barcode = in_barcode.split(';')
#     except AttributeError:  # means no fragments fall into the promoter
#         continue
#     tmp = pd.Series(in_barcode).value_counts().to_frame(name='count')
#     tmp['gene'] = np.full(shape=tmp.shape[0],fill_value=in_gene)
#     tmp.reset_index(inplace=True)  # three column: index, count, gene
#     bucket.append(tmp)
# df = pd.concat(bucket)
# df = df.loc[df['index'].isin(valid),:]
# final = df.groupby(by=['index','gene'])['count'].sum().unstack(fill_value=0)
# final.to_csv('atac_count_matrix.txt',sep='\t')

# count = pd.read_csv('atac_count_matrix.txt',sep='\t',index_col=0)
# def total_normalization(mat,target=1e4):
#     total = np.sum(mat,axis=1).reshape(-1,1)
#     sf = total/target
#     post = np.log(mat/sf + 1)
#     return post
# count = pd.DataFrame(data=total_normalization(count.values),index=count.index,columns=count.columns)
# count.to_csv('atac_count_norm.txt',sep='\t')

# count = pd.read_csv('atac_count_norm.txt',sep='\t',index_col=0)
# random = count.sample(n=30,axis=1)
# for gene in random.columns:
#     fig,ax = plt.subplots()
#     sns.histplot(random[gene].values,ax=ax)
#     ax.set_title('{}'.format(gene))
#     plt.savefig('dist/{}.pdf'.format(gene),bbox_inches='tight')
#     plt.close()


adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)
adata.var_names_make_unique()
adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
adata_peak = adata[:,adata.var['feature_types']=='Peaks']

adata_rna.write('pbmc_rna_to_azimuth.h5ad')
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
sc.pp.normalize_total(adata_rna,target_sum=1e4)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna,flavor='seurat',n_top_genes=3000)
adata_rna.raw = adata_rna
adata_rna = adata_rna[:,adata_rna.var['highly_variable']]
sc.pp.regress_out(adata_rna,['total_counts','pct_counts_mt'])
sc.pp.scale(adata_rna,max_value=10)
sc.tl.pca(adata_rna)
sc.pp.neighbors(adata_rna,n_pcs=40,n_neighbors=15)
sc.tl.leiden(adata_rna,resolution=1,key_added='rna_leiden1')
sc.tl.leiden(adata_rna,resolution=2,key_added='rna_leiden2')
sc.tl.leiden(adata_rna,resolution=3,key_added='rna_leiden3')
sc.tl.umap(adata_rna)

azimuth = pd.read_csv('./azimuth_pred.tsv',sep='\t',index_col=0)
azimuth_map = azimuth['predicted.celltype.l2'].to_dict()
azimuth_prediction = azimuth['predicted.celltype.l2.score'].to_dict()
azimuth_mapping = azimuth['mapping.score'].to_dict()
adata_rna.obs['azimuth'] = adata_rna.obs_names.map(azimuth_map).values
adata_rna.obs['prediction_score'] = adata_rna.obs_names.map(azimuth_prediction).values
adata_rna.obs['mapping_score'] = adata_rna.obs_names.map(azimuth_mapping).values

def final_plot(adata,col):
    fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(8,20),gridspec_kw={'hspace':0.3})  # for final_annotation
    sc.pl.umap(adata,color=col,frameon=False,ax=ax[0])
    sc.pl.umap(adata,color=col,frameon=False,legend_loc='on data',legend_fontsize=5,ax=ax[1])
    plt.savefig('umap_{}.pdf'.format(col),bbox_inches='tight')
    plt.close()

# final_plot(adata_rna,'rna_leiden1')
# final_plot(adata_rna,'rna_leiden2')
# final_plot(adata_rna,'rna_leiden3')
final_plot(adata_rna,'azimuth')

count = pd.read_csv('atac_count_norm.txt',sep='\t',index_col=0)
adata_active = ad.AnnData(X=count.values,var=pd.DataFrame(index=count.columns),obs=pd.DataFrame(index=count.index))
sc.pp.highly_variable_genes(adata_active,flavor='seurat',n_top_genes=3000)
adata_active.raw = adata_active
adata_active = adata_active[:,adata_active.var['highly_variable']]
sc.pp.scale(adata_active,max_value=10)
sc.tl.pca(adata_active)
sc.pp.neighbors(adata_active,n_pcs=40,n_neighbors=15)
sc.tl.leiden(adata_active,resolution=1,key_added='active_leiden1')
sc.tl.leiden(adata_active,resolution=2,key_added='active_leiden2')
sc.tl.leiden(adata_active,resolution=3,key_added='active_leiden3')
#sc.tl.umap(adata_active)
adata_active = adata_active[adata_rna.obs_names,:]   # make sure the order of cell is the same as rna
adata_active.obsm['X_umap'] = adata_rna.obsm['X_umap']

final_plot(adata_active,'active_leiden1')
final_plot(adata_active,'active_leiden2')
final_plot(adata_active,'active_leiden3')

# combine RNA and ATAC
adata_rna = adata_rna.raw.to_adata()
adata_active = adata_active.raw.to_adata()
assert np.all(adata_rna.obs_names.values == adata_active.obs_names.values)
print(adata_rna,adata_active)
print(adata_rna.X)
print(adata_active.X)
X = np.concatenate([adata_rna.X.toarray(),adata_active.X],axis=1)
var_names = adata_rna.var_names.to_list() + ['@@'+item for item in adata_active.var_names]
obs_names = adata_rna.obs_names
adata_combine = ad.AnnData(X=X,var=pd.DataFrame(index=var_names),obs=pd.DataFrame(index=obs_names.values))
adata_combine.obs['rna_leiden1'] = adata_rna.obs['rna_leiden1']
adata_combine.obs['rna_leiden2'] = adata_rna.obs['rna_leiden2']
adata_combine.obs['rna_leiden3'] = adata_rna.obs['rna_leiden3']
adata_combine.obs['active_leiden1'] = adata_active.obs['active_leiden1']
adata_combine.obs['active_leiden2'] = adata_active.obs['active_leiden2']
adata_combine.obs['active_leiden3'] = adata_active.obs['active_leiden3']
adata_combine.obsm['X_umap'] = adata_rna.obsm['X_umap']
print(adata_combine)
adata_combine.write('adata_combine.h5ad')









