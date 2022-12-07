#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os,sys
from sklearn.preprocessing import Binarizer
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import *


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# inspect data
adata_rna = sc.read_10x_h5('filtered_gene_bc_matrices_h5.h5')   # 10522 × 32738
adata_rna.var_names_make_unique()
# adata_full_splice = small_txt_to_adata('imputed-all-PSI-scaled.txt',True)   # 8665 × 17518
adata_full_splice = sc.read('adata_full_splice.h5ad')
adata_part_splice = small_txt_to_adata('imputed-all-PSI-scaled-400.txt',True)   # 8665 × 432
annotations = pd.read_csv('ICGS-splice-Subclusters.txt',sep='\t',index_col=0)


# run1: full_splice, all annotations
# common = list(set(adata_rna.obs_names).intersection(set(adata_full_splice.obs_names)).intersection(set(annotations.index)))  # 8665
# adata_rna = adata_rna[common,:]
# adata_rna.write('to_wnn_rna.h5ad')
# adata_rna = scanpy_recipe(adata_rna,False,[0.5,1,2,3,4,5],'rna',True,pca_n_comps=50,n_top_genes=3000)
# adata_rna = sc.read('full_adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_umap_True.h5ad')
#umap_dual_view_save(adata_rna,cols=['sctri_rna_leiden_{}'.format(r) for r in [0.5,1,2,3,4,5]])
# adata_full_splice = adata_full_splice[common,:]
# adata_full_splice.X = make_sure_mat_dense(adata_full_splice.X) * 10
# adata_full_splice.write('to_wnn_splice.h5ad')
# annotations = annotations.loc[common,:]
# adata_combined = concat_rna_and_other(adata_rna,adata_full_splice,'rna','splice','splice_')
# add_annotations(adata_combined,annotations,['ICGS2-Cluster','Final-Subcluster-Name'],kind='memory')
# # umap_dual_view_save(adata_combined,cols=['ICGS2-Cluster','Final-Subcluster-Name'])
# sctri = ScTriangulate(dir='output_full_two',adata=adata_combined,query=['ICGS2-Cluster','Final-Subcluster-Name'])
# sctri.lazy_run(viewer_heterogeneity_keys=['ICGS2-Cluster'],nca_embed=False)




# run2: part_splice, part_annotations
annotations = annotations.loc[[False if '--unassigned' in item else True for item in annotations['Final-Subcluster-Name']],:]
common = list(set(adata_rna.obs_names).intersection(set(adata_part_splice.obs_names)).intersection(set(annotations.index)))
adata_rna = adata_rna[common,:]  # 4013 × 32738
adata_part_splice = adata_part_splice[common,:]  # 4013 × 432
annotations = annotations.loc[common,:]
# # adata_rna = scanpy_recipe(adata_rna,False,[0.5,1,2,3,4,5],'rna',True,pca_n_comps=50,n_top_genes=3000)
adata_rna = sc.read('part_adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_umap_True.h5ad')
# # umap_dual_view_save(adata_rna,cols=['sctri_rna_leiden_{}'.format(r) for r in [0.5,1,2,3,4,5]])
adata_combined = concat_rna_and_other(adata_rna,adata_part_splice,'rna','splice','splice_')
add_annotations(adata_combined,annotations,['ICGS2-Cluster','Final-Subcluster-Name'],kind='memory')
# # umap_dual_view_save(adata_combined,cols=['ICGS2-Cluster','Final-Subcluster-Name'])
# sctri = ScTriangulate(dir='output_part_two',adata=adata_combined,query=['ICGS2-Cluster','Final-Subcluster-Name'])
# sctri.lazy_run(viewer_heterogeneity_keys=['ICGS2-Cluster'],nca_embed=False)
sctri = ScTriangulate.deserialize('output_part_two/after_pruned_assess.p')
# sc.pp.highly_variable_genes(adata_rna,flavor='seurat',n_top_genes=3000)
# hv_genes = adata_rna.var['highly_variable'].loc[adata_rna.var['highly_variable']].index.tolist()
# hv_splices = ['splice_' + item for item in adata_part_splice.var_names]
# hv_features = hv_genes + hv_splices
# adata_nca = nca_embedding(sctri.adata,10,'pruned','umap',hv_features=hv_features)
# adata_nca.write('output_part_two/adata_nca.h5ad')
adata_nca = sc.read('output_part_two/adata_nca.h5ad')
# ScTriangulate.salvage_run('build_all_viewers','output_part_two/after_pruned_assess.p',viewer_heterogeneity_keys=['ICGS2-Cluster'],other_umap=adata_nca.obsm['X_umap'])
# sctri.adata.obs['pruned'].to_csv('output_part_two/barcode2pruned.txt',sep='\t')
# pd.DataFrame(data=adata_nca.obsm['X_umap'],columns=['umap_x','umap_y'],index=adata_nca.obs_names).to_csv('output_part_two/nca_umap_coords.txt',sep='\t')
# sctri.modality_contributions(regex_dict={'splicing':r'^splice_'})
# for col in ['confidence','splicing_contribution','rna_contribution']:
#     sctri.adata.obsm['X_umap'] = adata_nca.obsm['X_umap']
#     sctri.plot_umap(col,'continuous',umap_cmap='viridis')
# for col in ['pruned','ICGS2-Cluster','final_annotation']:
#     sctri.plot_umap(col,'category')
# sctri.plot_multi_modal_feature_rank(cluster='Final-Subcluster-Name@MultiLin_c11--sub1',regex_dict={'splicing':r'^splice_'})
# sctri.plot_multi_modal_feature_rank(cluster='ICGS2-Cluster@MultiLin_c11',regex_dict={'splicing':r'^splice_'})
# sctri.plot_multi_modal_feature_rank(cluster='Final-Subcluster-Name@MKP_c4--sub3',regex_dict={'splicing':r'^splice_'})
sctri.adata.obsm['X_umap'] = adata_nca.obsm['X_umap']
# sc.pl.umap(sctri.adata,color=['splice_LAT2__ENSG00000086730__E2.10-E3.1&ENSG00000086730__E2.10-E4.1'],color_map=bg_greyed_cmap('viridis'),vmin=1e-5);plt.savefig('output_part_two/demo1_lat2.pdf',bbox_anchor='tight');plt.close()
# sctri.gene_to_df('marker_genes','pruned')
# ScTriangulate.salvage_run('build_all_viewers','output_part_two/after_pruned_assess.p',viewer_cluster=False,viewer_heterogeneity_keys=['ICGS2-Cluster'],other_umap=adata_nca.obsm['X_umap'],heatmap_scale=True,cmap=retrieve_pretty_cmap('altanalyze'))
# sctri.plot_two_column_sankey('ICGS2-Cluster','pruned',opacity=0.5,margin=300,text=False)
# sctri.plot_two_column_sankey('ICGS2-Cluster','pruned',opacity=0.5,margin=5,text=True)
# ScTriangulate.salvage_run('build_all_viewers','output_part_two/after_pruned_assess.p',viewer_cluster=False,viewer_heterogeneity_keys=['ICGS2-Cluster'],
#                            other_umap=adata_nca.obsm['X_umap'],heatmap_scale='median',heatmap_cmap=retrieve_pretty_cmap('altanalyze'),heatmap_regex=r'^splice_',
#                            heatmap_direction='exclude',heatmap_n_genes=8,heatmap_cbar_scale=0.25)
# ScTriangulate.salvage_run('build_all_viewers','output_part_two/after_pruned_assess.p',viewer_cluster=False,viewer_heterogeneity_keys=['ICGS2-Cluster'],
#                            other_umap=adata_nca.obsm['X_umap'],heatmap_scale='median',heatmap_cmap=retrieve_pretty_cmap('altanalyze'),heatmap_regex=r'^splice_',
#                            heatmap_direction='include',heatmap_n_genes=8,heatmap_cbar_scale=(-0.5,0.5))
# plot_coexpression(sctri.adata,'splice_LAT2__ENSG00000086730__E2.10-E3.1&ENSG00000086730__E2.10-E4.1',
#                   'splice_LAT2__ENSG00000086730__E2.10-E4.1&ENSG00000086730__E2.10-E3.1',kind='contour',contour_levels=20,outdir='output_part_two')
# sctri.plot_heterogeneity(key='ICGS2-Cluster',cluster='IER_c6',style='coexpression',gene1='splice_LAT2__ENSG00000086730__E2.10-E3.1&ENSG00000086730__E2.10-E4.1',
#                     gene2='splice_LAT2__ENSG00000086730__E2.10-E4.1&ENSG00000086730__E2.10-E3.1',kind='contour',contour_levels=20)
for cluster in ['MultiLin_c11','Cell_Cycle_c5','IER_c6','MKP_c4']:
    sctri.plot_heterogeneity(key='ICGS2-Cluster',cluster=cluster,style='build',heatmap_scale='median',heatmap_cmap=retrieve_pretty_cmap('altanalyze'),heatmap_regex=r'^splice_',heatmap_direction='exclude',heatmap_n_genes=8,heatmap_cbar_scale=0.25)
# wnn_meta = pd.read_csv('wnn_metadata.txt',sep='\t',index_col=0)
# add_annotations(sctri.adata,inputs=wnn_meta,kind='memory',cols_input=wnn_meta.columns.tolist())
# sctri.adata.obs['RNA.weight'] = sctri.adata.obs['RNA.weight'].astype(np.float64)
# sctri.adata.obs['SPLICE.weight'] = sctri.adata.obs['SPLICE.weight'].astype(np.float64)
# for col in ['RNA.weight','SPLICE.weight']:
#     sctri.plot_umap(col=col,kind='continuous',umap_cmap='viridis',umap_dot_size=4)
# sctri.cluster_performance(cluster=wnn_meta.columns.tolist()[-1],competitors=wnn_meta.columns.tolist()[2:-1],reference='pruned',show_cluster_number=True,ylim=(0,1))




