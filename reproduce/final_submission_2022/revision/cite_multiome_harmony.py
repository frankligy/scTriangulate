#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6


'''
old env: #!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6
new env: #!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7
old codebase: /data/salomonis2/LabFiles/Frank-Li/scTriangulate/src_backup
new codebase: /data/salomonis2/software
'''

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
import os,sys
sys.path.insert(0,'/data/salomonis2/LabFiles/Frank-Li/scTriangulate/src_backup')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap
from scipy.sparse import csr_matrix,issparse
import harmonypy as hm


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


'''
Demonstrating if we can run scTriangulate on a latent space from CITE and Multiome,
After having latent space, we first inspect how good it is, and then we can either

1. triangulate leiden resolutions
2. KNN transfer labels and then triangulate
'''


adata_cite = ScTriangulate.deserialize('after_pruned_assess_cite_tnc1.p').adata
adata_multiome = ScTriangulate.deserialize('after_pruned_assess_multiome.p').adata

'''
print(adata_cite)

AnnData object with n_obs × n_vars = 6406 × 32769
    obs: 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'n_genes', 'n_counts', 'doublet_scores', 'azimuth', 'prediction_score', 'mapping_score', 'sctri_rna_leiden_0.5', 'sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3', 'sctri_rna_leiden_4', 'sctri_rna_leiden_5', 'sctri_rna_leiden_6', 'n_adts_by_counts', 'sctri_adt_leiden_0.5', 'sctri_adt_leiden_1', 'sctri_adt_leiden_2', 'sctri_adt_leiden_3', 'sctri_adt_leiden_4', 'sctri_adt_leiden_5', 'sctri_adt_leiden_6', 'reassign@sctri_adt_leiden_1', 'tfidf10@sctri_adt_leiden_1', 'SCCAF@sctri_adt_leiden_1', 'doublet@sctri_adt_leiden_1', 'reassign@sctri_adt_leiden_2', 'tfidf10@sctri_adt_leiden_2', 'SCCAF@sctri_adt_leiden_2', 'doublet@sctri_adt_leiden_2', 'reassign@sctri_adt_leiden_3', 'tfidf10@sctri_adt_leiden_3', 'SCCAF@sctri_adt_leiden_3', 'doublet@sctri_adt_leiden_3', 'reassign@sctri_rna_leiden_1', 'tfidf10@sctri_rna_leiden_1', 'SCCAF@sctri_rna_leiden_1', 'doublet@sctri_rna_leiden_1', 'reassign@sctri_rna_leiden_2', 'tfidf10@sctri_rna_leiden_2', 'SCCAF@sctri_rna_leiden_2', 'doublet@sctri_rna_leiden_2', 'reassign@sctri_rna_leiden_3', 'tfidf10@sctri_rna_leiden_3', 'SCCAF@sctri_rna_leiden_3', 'doublet@sctri_rna_leiden_3', 'final_annotation', 'sctri_adt_leiden_1_shapley', 'sctri_adt_leiden_2_shapley', 'sctri_adt_leiden_3_shapley', 'sctri_rna_leiden_1_shapley', 'sctri_rna_leiden_2_shapley', 'sctri_rna_leiden_3_shapley', 'raw', 'prefixed', 'reassign@raw', 'tfidf10@raw', 'SCCAF@raw', 'doublet@raw', 'pruned', 'confidence', 'ori', 'reassign@pruned', 'tfidf10@pruned', 'SCCAF@pruned', 'doublet@pruned'
    var: 'gene_ids', 'feature_types', 'genome', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'modality'
    uns: 'azimuth_colors', 'sctri_adt_leiden_0.5_colors', 'sctri_adt_leiden_1_colors', 'sctri_adt_leiden_2_colors', 'sctri_adt_leiden_3_colors', 'sctri_adt_leiden_4_colors', 'sctri_adt_leiden_5_colors', 'sctri_adt_leiden_6_colors', 'sctri_rna_leiden_0.5_colors', 'sctri_rna_leiden_1_colors', 'sctri_rna_leiden_2_colors', 'sctri_rna_leiden_3_colors', 'sctri_rna_leiden_4_colors', 'sctri_rna_leiden_5_colors', 'sctri_rna_leiden_6_colors', 'final_annotation_colors', 'raw_colors', 'pruned_colors'
    obsm: 'X_pca', 'X_umap'
'''

'''
print(adata_multiome)

AnnData object with n_obs × n_vars = 10991 × 144945
    obs: 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'doublet_scores', 'azimuth', 'prediction_score', 'mapping_score', 'sctri_rna_leiden_0.5', 'sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3', 'sctri_rna_leiden_4', 'sctri_rna_leiden_5', 'sctri_rna_leiden_6', 'n_peaks_by_counts', 'sctri_atac_leiden_0.5', 'sctri_atac_leiden_1', 'sctri_atac_leiden_2', 'sctri_atac_leiden_3', 'sctri_atac_leiden_4', 'sctri_atac_leiden_5', 'sctri_atac_leiden_6', 'reassign@sctri_rna_leiden_1', 'tfidf10@sctri_rna_leiden_1', 'SCCAF@sctri_rna_leiden_1', 'doublet@sctri_rna_leiden_1', 'reassign@sctri_rna_leiden_2', 'tfidf10@sctri_rna_leiden_2', 'SCCAF@sctri_rna_leiden_2', 'doublet@sctri_rna_leiden_2', 'reassign@sctri_rna_leiden_3', 'tfidf10@sctri_rna_leiden_3', 'SCCAF@sctri_rna_leiden_3', 'doublet@sctri_rna_leiden_3', 'reassign@sctri_atac_leiden_1', 'tfidf10@sctri_atac_leiden_1', 'SCCAF@sctri_atac_leiden_1', 'doublet@sctri_atac_leiden_1', 'reassign@sctri_atac_leiden_2', 'tfidf10@sctri_atac_leiden_2', 'SCCAF@sctri_atac_leiden_2', 'doublet@sctri_atac_leiden_2', 'reassign@sctri_atac_leiden_3', 'tfidf10@sctri_atac_leiden_3', 'SCCAF@sctri_atac_leiden_3', 'doublet@sctri_atac_leiden_3', 'final_annotation', 'sctri_rna_leiden_1_shapley', 'sctri_rna_leiden_2_shapley', 'sctri_rna_leiden_3_shapley', 'sctri_atac_leiden_1_shapley', 'sctri_atac_leiden_2_shapley', 'sctri_atac_leiden_3_shapley', 'raw', 'prefixed', 'reassign@raw', 'tfidf10@raw', 'SCCAF@raw', 'doublet@raw', 'pruned', 'confidence', 'ori', 'reassign@pruned', 'tfidf10@pruned', 'SCCAF@pruned', 'doublet@pruned', 'tmp_plot', 'tfidf5@sctri_rna_leiden_1', 'tfidf5@sctri_rna_leiden_2', 'tfidf5@sctri_rna_leiden_3', 'tfidf5@sctri_atac_leiden_1', 'tfidf5@sctri_atac_leiden_2', 'tfidf5@sctri_atac_leiden_3', 'tfidf5@raw', 'tfidf5@pruned'
    var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'modality'
    uns: 'final_annotation_colors', 'raw_colors', 'pruned_colors'
    obsm: 'X_pca', 'X_umap'
'''

adata_rna_cite = adata_cite[:,adata_cite.var['feature_types']=='Gene Expression'].copy()
adata_rna_multiome = adata_multiome[:,adata_multiome.var['feature_types']=='Gene Expression'].copy()

common_genes = list(set(adata_rna_cite.var_names).intersection(set(adata_rna_multiome.var_names)))  # 20283
adata_rna_cite = adata_rna_cite[:,common_genes]
adata_rna_multiome = adata_rna_multiome[:,common_genes]

merge_rna_adata = ad.concat([adata_rna_cite,adata_rna_multiome],axis=0,join='outer',merge='first',label='batch',keys=['cite','multiome'])
merge_rna_adata.obs['concat'] = [b + ',' + a for b,a in zip(merge_rna_adata.obs['batch'],merge_rna_adata.obs['azimuth'])]
'''
print(merge_rna_adata)

AnnData object with n_obs × n_vars = 17397 × 20283
'''

make_sure_adata_writable(merge_rna_adata)
merge_rna_adata.write('merge_rna_adata.h5ad')


merge_rna_adata = sc.read('merge_rna_adata.h5ad')
def get_pca_for_rna(adata):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],percent_top=None,inplace=True,log1p=False)
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,flavor='seurat',n_top_genes=3000)
    adata.raw = adata
    adata = adata[:,adata.var['highly_variable']]
    sc.pp.regress_out(adata,['total_counts','pct_counts_mt'])
    sc.pp.scale(adata,max_value=10)
    sc.tl.pca(adata,n_comps=50)
    adata = adata.raw.to_adata()
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    return adata 
merge_rna_adata = get_pca_for_rna(merge_rna_adata)
ho = hm.run_harmony(merge_rna_adata.obsm['X_pca'],merge_rna_adata.obs,['batch'])
merge_rna_adata.obsm['X_pca_harmony'] = ho.Z_corr.T
sc.pp.neighbors(merge_rna_adata,use_rep='X_pca_harmony')
sc.tl.umap(merge_rna_adata)
umap_dual_view_save(merge_rna_adata,cols=['batch','azimuth','concat'])
merge_rna_adata.write('merge_rna_adata_harmony.h5ad')




