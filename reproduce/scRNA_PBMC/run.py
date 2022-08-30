#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import os
import sys
import scanpy as sc
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap,generate_gradient
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# this dataset is downloaded from
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3?
# and it is also available on synapse

# # load the data
# adata = sc.read_10x_h5('./pbmc_10k_v3_filtered_feature_bc_matrix.h5')
# adata.var_names_make_unique()
# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# adata.obs.to_csv('check_obs.txt',sep='\t')


# # qc
# for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
#     sc.pl.violin(adata,key,jitter=0.4)
#     plt.savefig('qc_violin_{}.pdf'.format(key),bbox_inches='tight')
#     plt.close()

# sc.pl.scatter(adata,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
# plt.savefig('qc_scatter.pdf',bbox_inches='tight')
# plt.close()


# # let's tentatively suggests 20% mitochondria, 500 min counts, 300 min genes
# print(adata)   # 11769 × 33538
# sc.pp.filter_cells(adata, min_genes=300)
# sc.pp.filter_cells(adata, min_counts=500)
# adata = adata[adata.obs.pct_counts_mt < 20, :]
# print(adata)   # 11022 × 33538
# adata.write('post_qc.h5ad')  # run Azimuth with this object

# # get the clustering and umap
# adata = sc.read('post_qc.h5ad')
# doublet_map = doublet_predict(adata)
# adding_azimuth(adata,'azimuth_pred.tsv')
# adata = scanpy_recipe(adata,is_log=False,resolutions=[0.5,1,2,3,4,5,6],pca_n_comps=50,n_top_genes=3000)
# adata.obs.to_csv('check_obs.txt',sep='\t')
# cols = adata.obs.columns.tolist()
# umap_dual_view_save(adata,cols)

# # run scTriangulate
adata = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
# adata.obs.to_csv('to_run_clustree.txt',sep='\t')
# sctri = ScTriangulate(dir='./output_two',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],predict_doublet='precomputed')
# sctri.lazy_run(viewer_heterogeneity_keys=[sctri.reference,'azimuth'])

# adata = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
# sctri = ScTriangulate(dir='./output_one',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],predict_doublet='precomputed',add_metrics={})
# sctri.lazy_run(viewer_heterogeneity_keys=[sctri.reference,'azimuth'])

# insights
# sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')

'''overall insights'''
# sctri.obs_to_df()

'''looking from leiden1'''
# # cluster0, CD4 TCM and CD4 TEM, sctriangulate success
# subset = ['sctri_rna_leiden_1@0','sctri_rna_leiden_3@2']
# sctri.plot_heterogeneity('sctri_rna_leiden_1','0','violin',subset=subset,genes=['CCR7','CCL5'])
# sctri.plot_heterogeneity('sctri_rna_leiden_1','0','dual_gene',subset=subset,dual_gene=['CCR7','CCL5'])
# sctri.plot_heterogeneity('sctri_rna_leiden_1','0','umap',subset=subset)

# cluster1,4,5,12,17, CD14 Monocyte, azimuth undercluster, sctriangulate identify new populations
# subset = ['sctri_rna_leiden_2@0','sctri_rna_leiden_3@10','sctri_rna_leiden_3@7','sctri_rna_leiden_1@4','sctri_rna_leiden_1@5','sctri_rna_leiden_3@13','sctri_rna_leiden_1@17']
# merge = [('sctri_rna_leiden_1@5','sctri_rna_leiden_3@13'),('sctri_rna_leiden_2@0','sctri_rna_leiden_3@10','sctri_rna_leiden_3@7','sctri_rna_leiden_1@4')]
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CD14',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='FCGR3A',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC10A',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC5A',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC4D',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='MX1',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='MX2',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='IFI44',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='S100A9',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='HLA-DRA',cmap='viridis')
# sctri.plot_heterogeneity('azimuth','CD14 Mono','umap',subset=subset)   
# sctri.plot_heterogeneity('azimuth','CD14 Mono','build',subset=subset)
# sctri.plot_heterogeneity('azimuth','CD14 Mono','umap',subset=subset,merge=merge)

# cluster2, CD4 naive, same

# cluster3, should be B memory + B intermediate, but none of the leiden recover it, I think it is due to the marker is not stable
# sc.pl.umap(sctri.adata,color=['MS4A1','TNFRSF13B','IGHM','IGHD','AIM2','CD79A','LINC01857','RALGPS2','BANK1','CD79B'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('check_b.pdf',bbox_inches='tight')
# plt.close()

# sc.pl.umap(sctri.adata,color=['MS4A1','COCH','AIM2','BANK1','SSPN','CD79A','TEX9','RALGPS2','TNFRSF13C','LINC01781'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('check_b1.pdf',bbox_inches='tight')
# plt.close()

# cluster6, MAIT + gammaT, success
# sctri.plot_heterogeneity('sctri_rna_leiden_1','6','umap')
# sctri.plot_heterogeneity('sctri_rna_leiden_1','6','dual_gene',dual_gene=['SLC4A10','TRDC'])

# cluster7, NK, cd56 fail
# sctri.plot_heterogeneity('sctri_rna_leiden_1','7','umap')
# sctri.plot_heterogeneity('sctri_rna_leiden_1','7','single_gene',single_gene='NCAM1',cmap='viridis')

# cluster8, CD8 TCM and TEM
# subset = ['sctri_rna_leiden_2@10','sctri_rna_leiden_3@18']
# sctri.plot_heterogeneity('sctri_rna_leiden_1','8','umap',subset=subset)
# sctri.plot_heterogeneity('sctri_rna_leiden_1','8','dual_gene',subset=subset,dual_gene=['IL7R','GZMH'])

# cluster9, B naive, same

# cluster10, CD16 Mono, same

# cluster11, CD8 naive

# cluster12, potentially doublet, predicted by scrublet with high confidence, also no obvious marker genes

# cluster13 has three parts, Treg, proliferating, plasmablast
# sctri.plot_heterogeneity('sctri_rna_leiden_1','13','umap')
# sctri.plot_heterogeneity('sctri_rna_leiden_1','13','multi_gene',multi_gene=['MZB1','IL2RA','MKI67'])

# cluster 14, potentially doublet, also no obvious marker genes

# cluster 15, platelet, same

# cluster 16, since we will remove ASDC, so that's fine, same, just cDC2

# cluster 17, potentially doublet, no obvious marker genes

# cluster 19, pDC, same

# cluster 20, ILC + HSPC
# sc.pl.umap(sctri.adata,color=['KIT','TRDC','TTLL10','LINC01229','SOX4','KLRB1','TNFRSF18','TNFRSF4', 'IL1R1', 'HPGDS'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('check_ilc.pdf',bbox_inches='tight')
# plt.close()

# cluster 21, cDC1, same


'''addtional plots for figure'''

# score_justify plot, using MAIT+gammaT as an example
# subset = ['sctri_rna_leiden_2@10','sctri_rna_leiden_3@18']
# sctri.plot_heterogeneity('sctri_rna_leiden_1','8','umap',subset=subset)
# subset = ['sctri_rna_leiden_3@18']
# sctri.plot_heterogeneity('sctri_rna_leiden_3','18','umap',subset=subset)
# sys.exit('stop')

# y1 = [0.69,0.41,0.86,0.49,0.02]
# y2 = [0.90,0.58,0.93,0.66,6.67]
# x1 = [(lambda x:3*x+1)(i) for i in range(5)]
# x2 = [(lambda x:3*x+2)(i) for i in range(5)]
# x = x1 + x2
# y = y1 + y2
# fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':[0.3,0.7],'hspace':0.1})
# for ax in axes:
#     ax.bar(x=x,height=y,color=np.repeat(['#1f77b4','#ff7f0e'],5),edgecolor='k')
#     for i in range(len(x)):
#         ax.text(x[i],y[i]+0.1,str(y[i]),va='center',ha='center')

# axes[0].set_ylim([6,7])
# axes[1].set_ylim([0,1])

# axes[0].spines['bottom'].set_visible(False)
# axes[1].spines['top'].set_visible(False)
# axes[0].tick_params(bottom=False)

# d = 0.015
# axes[0].plot((-d,d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
# axes[0].plot((1-d,1+d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
# axes[1].plot((-d,d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
# axes[1].plot((1-d,1+d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')

# t = [(lambda x:3*x+1.5)(i) for i in range(5)]
# axes[1].set_xticks(t)
# axes[1].set_xticklabels(['Reassign Score','SCCAF Score','TF-IDF10 Score','TF-IDF 5 Score','Shapley Value'],rotation=60)

# axes[1].legend(handles=[mpatch.Patch(color=i) for i in ['#1f77b4','#ff7f0e']],labels=['leiden1@9','leiden3@7'],
#           loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# axes[1].set_xlabel('Metrics and Shapley')
# axes[1].set_ylabel('Value')

# plt.savefig('score_justify.pdf',bbox_inches='tight')
# plt.close()

# confidence plot
# sctri.plot_umap('confidence','continuous',umap_cmap='viridis')

# pruned_final_annotation
# sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
# sctri.plot_umap('pruned_final_annotation','category')

# plot concordance
# confusion_df = sctri.plot_concordance(key1='azimuth',key2='pruned',style='3dbar')
# confusion_df.to_csv('output_two/azimuth_pruned_confusion.txt',sep='\t')



# cluster performance plot
# sctri.cluster_performance(cluster='pruned',competitors=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],reference='azimuth',show_cluster_number=True,metrics=None)



'''for building the tutorial'''
# sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
# subset = ['sctri_rna_leiden_2@7','sctri_rna_leiden_2@18']
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='umap',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='build',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='violin',genes=['TRDC','SLC4A10'],subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='single_gene',single_gene='TRDC',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='dual_gene',dual_gene=['TRDC','SLC4A10'],subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='sankey',subset=subset)

# sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1')
# sctri.plot_circular_barplot(key='sctri_rna_leiden_1',col='pruned')

# sctri.get_metrics_and_shapley(barcode='AAACCCACATCCAATG-1',save=True)

# generate_gradient(cmap='jet',name='jet')
# generate_block(color_list=pick_n_colors(200),name='r433')


'''
revision sensitivity test
'''

# for opt in ['rank_all_or_none','rank','shapley']:
#     ScTriangulate.salvage_run(step_to_start='run_shapley',last_step_file='output_two/after_metrics.p',outdir='output_two_{}'.format(opt),
#                               shapley_mode=opt,shapley_bonus=0.01,assess_pruned=False,viewer_cluster=False,viewer_heterogeneity=False)

# sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
# for opt in ['rank_all_or_none','rank','shapley']:
#     new_sctri = ScTriangulate.deserialize('output_two_{}/after_rank_pruning.p'.format(opt))
#     new_sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(new_sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
#     new_sctri.add_to_invalid_by_win_fraction(percent=0.25)
#     new_sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=new_sctri.reference)
#     new_obs = new_sctri.adata.obs
#     add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_{}'.format(opt)],'\t','memory')

# sctri.cluster_performance(cluster='pruned',competitors=['pruned_{}'.format(opt) for opt in ['rank_all_or_none','rank','shapley']],
#                           reference='azimuth',show_cluster_number=True,metrics=True)



# for bonus in [0,0.005,0.01,0.05,0.1]:
#     ScTriangulate.salvage_run(step_to_start='run_shapley',last_step_file='output_two/after_metrics.p',outdir='output_two_bonus_{}'.format(bonus),
#                               shapley_mode='shapley_all_or_none',shapley_bonus=bonus,assess_pruned=False,viewer_cluster=False,viewer_heterogeneity=False)

# sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
# for bonus in [0,0.005,0.01,0.05,0.1]:
#     new_sctri = ScTriangulate.deserialize('output_two_bonus_{}/after_rank_pruning.p'.format(bonus))
#     new_sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(new_sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
#     new_sctri.add_to_invalid_by_win_fraction(percent=0.25)
#     new_sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=new_sctri.reference)
#     new_obs = new_sctri.adata.obs
#     add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_bonus_{}'.format(bonus)],'\t','memory')

# sctri.cluster_performance(cluster='pruned',competitors=['pruned_bonus_{}'.format(bonus) for bonus in [0,0.005,0.01,0.05,0.1]],
#                           reference='azimuth',show_cluster_number=True,metrics=True)


# sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
# for wfc in [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]:
#     sctri.clear_invalid()
#     sctri.add_to_invalid_by_win_fraction(percent=wfc)
#     sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=sctri.reference)
#     new_obs = sctri.adata.obs
#     add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_wfc_{}'.format(wfc)],'\t','memory')

# sctri.adata.obs['pruned'] = sctri.adata.obs['pruned_wfc_0.25']
# sctri.cluster_performance(cluster='pruned',competitors=['pruned_wfc_{}'.format(wfc) for wfc in [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]],
#                           reference='azimuth',show_cluster_number=True,metrics=True)


sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
for nabs in [1,5,10,15,20,25,30,35,40,45,50]:
    sctri.clear_invalid()
    sctri.add_to_invalid_by_win_fraction(percent=0.25)
    sctri.pruning(method='reassign',abs_thresh=nabs,remove1=True,reference=sctri.reference)
    new_obs = sctri.adata.obs
    add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_nabs_{}'.format(nabs)],'\t','memory')

sctri.adata.obs['pruned'] = sctri.adata.obs['pruned_nabs_10']
sctri.cluster_performance(cluster='pruned',competitors=['pruned_nabs_{}'.format(nabs) for nabs in [1,5,10,15,20,25,30,35,40,45,50]],
                          reference='azimuth',show_cluster_number=True,metrics=True)






















