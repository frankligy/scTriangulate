#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import os
import sys
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import re
from sctriangulate import ScTriangulate
from sctriangulate.colors import bg_greyed_cmap
from sctriangulate.preprocessing import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

'''first part of the analysis is just relying on count matrix, no fragments involved
   to find closest genes, the gtf annotation is needed'''

# adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)  # 11909
# adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
# adata_atac = adata[:,adata.var['feature_types']=='Peaks']
# adata_atac = reformat_peak(adata_atac)
# find_genes(adata_atac,'gencode.v38.annotation.gtf')
# adata_atac.var_names = [str(v) + '_' + adata_atac.var['gene_annotation'][i] for i,v in enumerate(adata_atac.var_names)]


# QC, RNA
# adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
# sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
#     sc.pl.violin(adata_rna,key,jitter=0.4)
#     plt.savefig('qc_rna_violin_{}.pdf'.format(key),bbox_inches='tight')
#     plt.close()

# sc.pl.scatter(adata_rna,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
# plt.savefig('qc_rna_scatter.pdf',bbox_inches='tight')
# plt.close()

# QC, ATAC
# sc.pp.calculate_qc_metrics(adata_atac, var_type='peaks',percent_top=None, log1p=False, inplace=True)
# for key in ['n_peaks_by_counts','total_counts']:
#     sc.pl.violin(adata_atac,key,jitter=0.4)
#     plt.savefig('qc_atac_violin_{}.pdf'.format(key),bbox_inches='tight')
#     plt.close()

# sc.pl.scatter(adata_atac,x='n_peaks_by_counts',y='total_counts')
# plt.savefig('qc_peak_scatter.pdf',bbox_inches='tight')
# plt.close()

# decision
# for RNA, min_gene = 300, min_counts = 500, mito < 20%
# for ATAC, min_peak = 1000
# adata_rna = adata_rna[adata_rna.obs['n_genes_by_counts']>300,:]
# adata_rna = adata_rna[adata_rna.obs['total_counts']>500,:]
# adata_rna = adata_rna[adata_rna.obs['pct_counts_mt']<20,:]
# kept_rna = adata_rna.obs_names
# adata_atac = adata_atac[adata_atac.obs['n_peaks_by_counts']>1000,:]
# kept_atac = adata_atac.obs_names
# common = list(set(kept_rna).intersection(set(kept_atac)))

# adata_rna = adata_rna[adata_rna.obs_names.isin(common),:]   # 10991 × 36601
# adata_atac = adata_atac[adata_atac.obs_names.isin(common),:]  # 10991 × 108344
# adata_rna.write('post_rna_qc.h5ad')
# adata_atac.write('post_atac_qc.h5ad')


# get the clusters and combine
# doublet_map = doublet_predict(adata_rna)
# adding_azimuth(adata_rna,'azimuth_pred.tsv')
# adata_rna = scanpy_recipe(adata_rna,False,resolutions=[0.5,1,2,3,4,5,6],modality='rna',pca_n_comps=50,n_top_genes=3000)
# adata_atac = scanpy_recipe(adata_atac,False,resolutions=[0.5,1,2,3,4,5,6],modality='atac',pca_n_comps=100,n_top_genes=100000)

adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_atac = sc.read('adata_after_scanpy_recipe_atac_0.5_1_2_3_4_5_6_umap_True.h5ad')
# adata_combine = concat_rna_and_other(adata_rna,adata_atac,umap='other',name='atac',prefix='')
# adata_combine = make_sure_adata_writable(adata_combine,delete=True)
# adata_combine.write('combined_rna_atac_au.h5ad')


# all_r = ['sctri_{}_leiden_{}'.format(m,i) for i in [0.5,1,2,3,4,5,6] for m in ['rna','atac']]
# cols = all_r + ['azimuth','doublet_scores']
# umap_dual_view_save(adata_combine,cols)


# run scTriangulate
adata_combine = sc.read('combined_rna_atac_ru.h5ad')
# sctri = ScTriangulate(dir='output_one',adata=adata_combine,add_metrics={},query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3',
#                       'sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'],predict_doublet='precomputed')
# sctri.lazy_run(scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_atac_leiden_1'])

# sctri = ScTriangulate(dir='output_two',adata=adata_combine,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3',
#                       'sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'],predict_doublet='precomputed')
# sctri.lazy_run(scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_atac_leiden_1'])

'''I found out one tfidf is not enough, have to use two, but 25% cutoff for two is too messy,
   so, I first check which cutoff maybe better, I tried 25% - 40%, I found 35% seems reasonable,
   now I need to conduct salvage run'''
# ScTriangulate.salvage_run(step_to_start='assess_pruned',last_step_file='output_two/after_rank_pruning.p',win_fraction_cutoff=0.35)


# #### Re-run
# sctri = ScTriangulate(dir='output_two_new',adata=adata_combine,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3',
#                       'sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'],predict_doublet='precomputed')
# sctri.lazy_run(compute_metrics_parallel=False,scale_sccaf=True,win_fraction_cutoff=0.35,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_atac_leiden_1'])
# sys.exit('stop')


## insights
sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
# sctri.cluster_performance(cluster='pruned',competitors=sctri.query,reference='azimuth',show_cluster_number=True)
# cols = ['wsnn_res.0.25','wsnn_res.0.5','wsnn_res.0.75','wsnn_res.1','wsnn_res.1.25','wsnn_res.1.5','wsnn_res.1.75',
#         'wsnn_res.2','wsnn_res.2.25','wsnn_res.2.5','wsnn_res.2.75','wsnn_res.3','wsnn_res.3.25','wsnn_res.3.5','wsnn_res.3.75']
# add_annotations(sctri.adata,'./wnn_metadata.txt',cols,0,cols)
# sctri.cluster_performance(cluster='pruned',competitors=cols,reference='azimuth',show_cluster_number=True)

# sctri_atac_leiden_2@23
# sctri.plot_long_heatmap(figsize=(10,12))
# sctri.plot_umap('AFF3','continuous',umap_cmap='viridis')
sctri.plot_umap('CD14','continuous',umap_cmap='viridis')
sys.exit('stop')


# sctri.modality_contributions()
# for col in ['rna_contribution','atac_contribution']:
#    sctri.plot_umap(col,'continuous',umap_cmap='viridis')
# sctri.plot_umap('confidence','continuous',umap_cmap='viridis')


## CD4 populations
# for cluster in ['sctri_rna_leiden_3@5','sctri_rna_leiden_2@14','sctri_rna_leiden_1@13','sctri_rna_leiden_3@24','sctri_rna_leiden_3@23','sctri_rna_leiden_3@1','sctri_rna_leiden_2@21']:
#    sctri.plot_multi_modal_feature_rank(cluster=cluster)

## let's focus on the CD4 population and make some figures
# sctri.adata.obs['dummy_key'] = np.full(shape=sctri.adata.obs.shape[0],fill_value='dummy_cluster')
# subset = ['sctri_rna_leiden_2@21','sctri_rna_leiden_3@1','sctri_rna_leiden_3@23','sctri_rna_leiden_3@5','sctri_rna_leiden_3@24','sctri_rna_leiden_2@14','sctri_rna_leiden_1@13']
# sctri.plot_heterogeneity('dummy_key','dummy_cluster','umap','pruned',subset=subset)
# for gene in ['PTPN13','SLC4A10','IL2RA','ITGB1','FHIT','TSHZ2','chr20_52971902_52974673_TSHZ2','chr6_22045235_22045691_CASC15','chr6_22146145_22147503_NBAT1;CASC15','SOX4','chr6_22148300_22148794_NBAT1;CASC15','chr6_21876998_21877505_CASC15','CASC15']:
#    sctri.plot_heterogeneity('dummy_key','dummy_cluster','single_gene','pruned',subset=subset,single_gene=gene,cmap='viridis')
# sys.exit('stop')
# subset = ['CD4 Naive','CD4 TCM','Treg','CD4 TEM','MAIT']
# sctri.plot_heterogeneity('dummy_key','dummy_cluster','umap','azimuth',subset=subset)


# sctri.plot_umap('CASC15','continuous',umap_cmap='viridis')
# sctri.plot_umap('MAP3K4','continuous',umap_cmap='viridis')
# genes = ['CXCR6','GZMB','PRDM1','CXCR5','IL21','BCL6','CCR7','IL7R','S1PR1']
# sc.pl.umap(sctri.adata,color=genes,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('cd4_markers.pdf',bbox_inches='tight')
# plt.close()

# CD14
for gene in ['HLA-DRA','S100A9','FCGR3A','CLEC4D','CLEC5A','CLEC10A','VCAN','JAK2','DYSF']:
   sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',single_gene=gene,cmap='viridis')



































## insights
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
# sctri.modality_contributions()
# for col in ['atac_contribution','rna_contribution']:
#    sctri.plot_umap(col,'continuous')

# sctri.gene_to_df(mode='marker_genes',key='pruned')
# sctri.gene_to_df(mode='exclusive_genes',key='pruned')

# sc.pl.umap(sctri.adata,color=['CLEC10A','CLEC5A','CLEC4D','FCGR3A','S100A9','HLA-DRA'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# sc.pl.umap(sctri.adata,color=['chr5_756911_759750_ZDHHC11B','ZDHHC11B','chr6_161076038_161077022_MAP3K4','MAP3K4','chr14_64755519_64756581_SPTB','SPTB','chr7_70787751_70788862_AUTS2','AUTS2','chr10_33114209_33116813_IATPR','IATPR','chr21_45223468_45225379_ADARB1','ADARB1','chr1_57307181_57308005_DAB1','DAB1'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('check_features.pdf',bbox_inches='tight')
# plt.close()

# looking at mono
# for cluster in ['sctri_atac_leiden_1@0','sctri_rna_leiden_2@12','sctri_rna_leiden_1@7','sctri_atac_leiden_1@13','sctri_rna_leiden_3@13','sctri_rna_leiden_3@12','sctri_rna_leiden_2@15','sctri_rna_leiden_3@35','sctri_rna_leiden_1@8','sctri_atac_leiden_2@13']:
#    sctri.plot_multi_modal_feature_rank(cluster=cluster)

# looking at B
# for cluster in ['sctri_atac_leiden_1@6','sctri_atac_leiden_1@14','sctri_atac_leiden_2@23','sctri_atac_leiden_2@16']:
#     sctri.plot_multi_modal_feature_rank(cluster=cluster)

# looking at CD8 and CD4
# for cluster in ['sctri_atac_leiden_1@1','sctri_rna_leiden_1@13','sctri_rna_leiden_2@1','sctri_rna_leiden_2@0','sctri_atac_leiden_3@4','sctri_rna_leiden_2@11']:
#    sctri.plot_multi_modal_feature_rank(cluster=cluster)

'''LOOK AT AZIMUTH population'''
# CD16 mono
sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
sctri.plot_heterogeneity('azimuth','CD16 Mono','heatmap')
sctri.plot_heterogeneity('azimuth','CD16 Mono','umap')

sys.exit('stop')



'''
   second part is to derive gene activity score, now we need fragments file
   1. please refer to gene_activity_count_matrix_old_10x function docstring to see how to derive fall_in_promoter and fall_in_gene,
      but actually run gene_activity_count_matrix_new_10x, new cellranger pipeline allows multiple count for one fragment
   2. the atac_count.txt is what we get from here
'''

# mainly about how to compare the gene activity score versus gene
atac_count = pd.read_csv('atac_count.txt',sep='\t',index_col=0)
atac_count_norm = pd.DataFrame(data=Normalization.total_normalization(atac_count.values),index=atac_count.index,columns=atac_count.columns)
adata_activity = ad.AnnData(X=atac_count_norm.values,var=pd.DataFrame(index=atac_count_norm.columns),obs=pd.DataFrame(index=atac_count_norm.index))
adata_activity.X = make_sure_mat_sparse(adata_activity.X)
adata_activity.var_names = ['@@'+item for item in adata_activity.var_names]

features = ['MS4A1','@MS4A1','CD79A','@@CD79A','NKG7','@@NKG7','CD8A','@@CD8A']
    







# adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)
# atac_count = gene_activity_count_matrix_new_10x('fall_in_promoter.bed','fall_in_gene.bed',adata.obs_names.tolist())
# atac_count.to_csv('atac_count.txt',sep='\t')

atac_count = pd.read_csv('atac_count.txt',sep='\t',index_col=0)
atac_count_norm = pd.DataFrame(data=Normalization.total_normalization(atac_count.values),index=atac_count.index,columns=atac_count.columns)
adata_atac = ad.AnnData(X=atac_count_norm.values,var=pd.DataFrame(index=atac_count_norm.columns),obs=pd.DataFrame(index=atac_count_norm.index))
adata_atac.X = make_sure_mat_sparse(adata_atac.X)
adata_atac = scanpy_recipe(adata_atac,is_log=True,resolutions=[1,2,3],modality='atac',umap=True,save=True,pca_n_comps=100)
adata_atac = sc.read('adata_after_scanpy_recipe_atac_1_2_3_umap_True.h5ad')
adding_azimuth(adata_atac,result='./azimuth_pred.tsv')
umap_dual_view_save(adata_atac,cols=['azimuth','sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'])

sys.exit('stop')

# adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)
# adata.var_names_make_unique()
# adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
# adding_azimuth(adata_rna,result='./azimuth_pred.tsv')
# adata_rna = scanpy_recipe(adata_rna,is_log=False,resolutions=[1,2,3],modality='rna',umap=True,save=True)
# print(adata_rna)

# adata_atac = sc.read('adata_after_scanpy_recipe_atac_1_2_3_umap_True.h5ad')
# adata_rna = sc.read('adata_after_scanpy_recipe_rna_1_2_3_umap_True.h5ad')

# adata = concat_rna_and_other(adata_rna,adata_atac,umap=None,name='atac',prefix='@@')
# adata.obsm['X_umap_rna'] = adata_rna.obsm['X_umap']
# adata.obsm['X_umap_atac'] = adata_atac.obsm['X_umap']
# adata.obsm['X_umap'] = adata.obsm['X_umap_rna']
# print(adata)
# print(adata.var_names)
# adata.write('atac_rna_combine.h5ad')

# umap_dual_view_save(adata,cols=['azimuth','sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3','sctri_atac_leiden_1',
#                                 'sctri_atac_leiden_2','sctri_atac_leiden_3'])


# # sanity check of my ATAC pipeline
# adata = sc.read('atac_rna_combine.h5ad')
# sc.pl.umap(adata,color=['MS4A1','@@MS4A1','CD79A','@@CD79A','NKG7','@@NKG7','CD8A','@@CD8A'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('sanity_check.pdf',bbox_inches='tight')


# run sctriangulate only rna
# adata = sc.read('adata_after_scanpy_recipe_rna_1_2_3_umap_True.h5ad')
# sctri = ScTriangulate(dir='./output_rna',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],
#                         reference='sctri_rna_leiden_1',species='human',criterion=2,verbose=2)
# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})

# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=True)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley()
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_rank_pruning.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_rna/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()

# sctri = ScTriangulate.deserialize('output_rna/after_rank_pruning.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.25)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# # for col in ['final_annotation','raw','pruned']:
# #     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='azimuth')
# sctri.viewer_heterogeneity_html(key='azimuth')
# sys.exit('stop')

# run sctriangulate only rna + atac1
# adata = sc.read('atac_rna_combine.h5ad')
# sctri = ScTriangulate(dir='./output_rna_atac1',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3','sctri_atac_leiden_1'],
#                         reference='sctri_rna_leiden_1',species='human',criterion=2,verbose=2)
# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})

# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=True)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley()
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_rank_pruning.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_rna_atac1/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()

# sctri = ScTriangulate.deserialize('output_rna_atac1/after_rank_pruning.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.25)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# for col in ['final_annotation','raw','pruned']:
#     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='azimuth')
# sctri.viewer_heterogeneity_html(key='azimuth')
# sys.exit('stop')


# run sctriangulate only rna + atac12
# adata = sc.read('atac_rna_combine.h5ad')
# sctri = ScTriangulate(dir='./output_rna_atac12',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3','sctri_atac_leiden_1','sctri_atac_leiden_2'],
#                         reference='sctri_rna_leiden_1',species='human',criterion=2,verbose=2)
# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})

# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=True)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley()
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_rank_pruning.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_rna_atac12/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()

# sctri = ScTriangulate.deserialize('output_rna_atac12/after_rank_pruning.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.25)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# for col in ['final_annotation','raw','pruned']:
#     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='azimuth')
# sctri.viewer_heterogeneity_html(key='azimuth')
# sys.exit('stop')

#run sctriangulate only rna + atac123
# adata = sc.read('atac_rna_combine.h5ad')
# sctri = ScTriangulate(dir='./output_rna_atac123',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3','sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'],
#                         reference='sctri_rna_leiden_1',species='human',criterion=2,verbose=2)
# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})

# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=True)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley()
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_rank_pruning.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_rna_atac123/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()

# sctri = ScTriangulate.deserialize('output_rna_atac123/after_rank_pruning.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.25)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# for col in ['final_annotation','raw','pruned']:
#     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='azimuth')
# sctri.viewer_heterogeneity_html(key='azimuth')


# insight check
sctri_atac = ScTriangulate.deserialize('output_rna_atac1/after_rank_pruning.p')
sctri_atac.plot_heterogeneity('azimuth','Eryth','pruned','violin',save,'pdf',['@@PTPRD','PTPRD'])
sys.exit()

sctri_rna = ScTriangulate.deserialize('output_rna/after_rank_pruning.p')
mat = sctri_rna.adata.X.toarray()
print(mat.sum())
print(np.count_nonzero(mat)/mat.size)
# sctri_rna.add_to_invalid_by_win_fraction(percent=0.25)
# sctri_rna.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# cd4_tcm_marker_dict = sctri_rna.plot_heterogeneity('azimuth','CD4 TCM','pruned','heatmap',False)
# for cluster in ['sctri_rna_leiden_2@1','sctri_rna_leiden_3@8','sctri_rna_leiden_1@3','sctri_rna_leiden_3@18']:
#     cd4_tcm_marker_dict[cluster].to_csv('output_rna/{}_marker.txt'.format(cluster),sep='\t')

sctri_atac = ScTriangulate.deserialize('output_rna_atac12/after_rank_pruning.p')
mat = sctri_atac.adata.X.toarray()
print(mat.sum())
print(np.count_nonzero(mat)/mat.size)
# sc.pl.umap(sctri_atac.adata,color=['MEIS1','@@MEIS1'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('atac_l1_c9.pdf',bbox_inches='tight')
# sc.pl.umap(sctri_atac.adata,color=['DCC','@@DCC','TMEM132D','@@TMEM132D','UNC5D','@@UNC5D','DPP6','@@DPP6','KCNH5','@@KCNH5'],
#            cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('atac_l1_c7.pdf',bbox_inches='tight')
#sc.pl.umap(sctri_atac.adata,color=['KRT73','@@KRT73','ISM1','@@ISM1','EDAR','@@EDAR','WNT7A','@@WNT7A'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('cd8_naive_cd4_naive_cd4_tcm_umap.pdf',bbox_inches='tight')
# sctri_atac.add_to_invalid_by_win_fraction(percent=0.25)
# sctri_atac.pruning(method='reassign',abs_thresh=10,remove1=True,reference='azimuth')
# cd4_tcm_marker_dict = sctri_atac.plot_heterogeneity('azimuth','CD4 TCM','pruned','heatmap',False)
# for cluster in ['sctri_rna_leiden_2@1','sctri_rna_leiden_3@8','sctri_rna_leiden_1@3','sctri_rna_leiden_3@18']:
#     cd4_tcm_marker_dict[cluster].to_csv('output_rna_atac12/{}_marker.txt'.format(cluster),sep='\t')




























