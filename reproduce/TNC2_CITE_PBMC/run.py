#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import pandas as pd
import numpy as np
import os,sys
import scanpy as sc
import matplotlib.pyplot as plt

from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap

# load the data
# adata = sc.read_10x_h5('31WF_ND19-446__TNC-RNA-ADT.h5',gex_only=False)
# adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
# adata_adt = adata[:,adata.var['feature_types']=='Antibody Capture']

# adata_rna.var_names_make_unique()
# adata_adt.var_names_make_unique()


# qc, rna
# adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
# sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
#     sc.pl.violin(adata_rna,key,jitter=0.4)
#     plt.savefig('qc_rna_violin_{}.pdf'.format(key),bbox_inches='tight')
#     plt.close()

# sc.pl.scatter(adata_rna,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
# plt.savefig('qc_rna_scatter.pdf',bbox_inches='tight')
# plt.close()

# # qc, adt
# sc.pp.calculate_qc_metrics(adata_adt, var_type='adts',percent_top=None, log1p=False, inplace=True)

# for key in ['n_adts_by_counts','total_counts']:
#     sc.pl.violin(adata_adt,key,jitter=0.4)
#     plt.savefig('qc_adt_violin_{}.pdf'.format(key),bbox_inches='tight')
#     plt.close()

# sc.pl.scatter(adata_adt,x='n_adts_by_counts',y='total_counts')
# plt.savefig('qc_adt_scatter.pdf',bbox_inches='tight')
# plt.close()

# decision, min_genes = 300, min_counts= 500, mito < 20
# sc.pp.filter_cells(adata_rna, min_genes=300)
# sc.pp.filter_cells(adata_rna, min_counts=500)
# adata_rna = adata_rna[adata_rna.obs.pct_counts_mt < 20, :]
# adata_adt = adata_adt[adata_rna.obs_names,:]

# get clusters and input
# doublet_map = doublet_predict(adata_rna)
# adding_azimuth(adata_rna,'azimuth_pred.tsv')
# adata_rna = scanpy_recipe(adata_rna,False,resolutions=[0.5,1,2,3,4,5,6],modality='rna',pca_n_comps=50)
# adata_adt = scanpy_recipe(adata_adt,False,resolutions=[0.5,1,2,3,4,5,6],modality='adt',pca_n_comps=15)

adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_adt = sc.read('adata_after_scanpy_recipe_adt_0.5_1_2_3_4_5_6_umap_True.h5ad')
# adata_combine = concat_rna_and_other(adata_rna,adata_adt,umap='other',name='adt',prefix='AB_')

# all_r = ['sctri_{}_leiden_{}'.format(m,i) for i in [0.5,1,2,3,4,5,6] for m in ['rna','adt']]
# cols = all_r + ['azimuth','doublet_scores']
# umap_dual_view_save(adata_combine,cols)

# adata_combine = make_sure_adata_writable(adata_combine,delete=True)
# adata_combine.write('combined_rna_adt.h5ad')


# run scTriangulate
adata_combine = sc.read('combined_rna_adt.h5ad')
# sctri = ScTriangulate(dir='output_one',adata=adata_combine,add_metrics={},predict_doublet='precomputed',
#                       query=['sctri_adt_leiden_1','sctri_adt_leiden_2','sctri_adt_leiden_3','sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
# sctri.lazy_run(viewer_heterogeneity_keys=['azimuth','sctri_adt_leiden_1','sctri_rna_leiden_1'])


# sctri = ScTriangulate(dir='output_two',adata=adata_combine,predict_doublet='precomputed',
#                       query=['sctri_adt_leiden_1','sctri_adt_leiden_2','sctri_adt_leiden_3','sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
# sctri.lazy_run(viewer_heterogeneity_keys=['azimuth','sctri_adt_leiden_1','sctri_rna_leiden_1'])



# insights
## plot all adts
# all_adts = ['AB_' + item for item in adata_adt.var_names.tolist()]
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
# sc.pl.umap(sctri.adata,color=all_adts,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('all_adt.pdf',bbox_inches='tight')
# plt.close()

# sctri.plot_umap('AB_CD56','continuous',umap_cmap='viridis')

# cluster performance
# sctri.cluster_performance(cluster='pruned',competitors=sctri.query,reference='azimuth',show_cluster_number=True)
cols = ['wsnn_res.0.25','wsnn_res.0.5','wsnn_res.0.75','wsnn_res.1','wsnn_res.1.25','wsnn_res.1.5','wsnn_res.1.75',
        'wsnn_res.2','wsnn_res.2.25','wsnn_res.2.5','wsnn_res.2.75','wsnn_res.3','wsnn_res.3.25','wsnn_res.3.5','wsnn_res.3.75']
add_annotations(sctri.adata,'./wnn_metadata.txt',cols,0,cols)
sctri.cluster_performance(cluster='pruned',competitors=cols,reference='azimuth',show_cluster_number=True)

## output the obs
# sctri.obs_to_df()
# sys.exit('stop')

## adt and rna contributions
# sctri.modality_contributions()
# for col in ['adt_contribution','rna_contribution']:
#     sctri.plot_umap(col,'continuous',umap_cmap='viridis')


## distribution of resolutions
# col = []
# for item in sctri.adata.obs['pruned']:
#     if 'leiden_1@' in item:
#         col.append('resolution1')
#     elif 'leiden_2@' in item:
#         col.append('resolution2')
#     elif 'leiden_3@' in item:
#         col.append('resolution3')
# sctri.adata.obs['resolution_distribution'] = col
# sctri.plot_umap('resolution_distribution','category')

## WNN
# wnn = pd.read_csv('wnn_metadata.txt',sep='\t',index_col=0)
# for col in ['wsnn_res.0.25','wsnn_res.1','wsnn_res.2','wsnn_res.3','wsnn_res.3.75']:
#     mapping = wnn[col].to_dict()
#     sctri.adata.obs[col] = [str(item) for item in sctri.adata.obs_names.map(mapping).values]
#     sctri.plot_umap(col,'category')

## change to WNN umap, see if the CD4 population looks cleaner
# add_umap(sctri.adata,'wnn_umap.txt',mode='pandas',cols=['wnnuMAP_1','wnnuMAP_2'])
# sctri.plot_umap('pruned','category')


## CD4 memory, plot rank
# for cluster in ['sctri_adt_leiden_2@10','sctri_adt_leiden_1@12','sctri_adt_leiden_2@22','sctri_adt_leiden_3@1','sctri_rna_leiden_2@2','sctri_rna_leiden_1@0','sctri_adt_leiden_1@0']:
#     sctri.plot_multi_modal_feature_rank(cluster=cluster)



