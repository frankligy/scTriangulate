#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os
import sys
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import *

'''
all the input files, intermediate results are deposited on https://www.synapse.org/#!Synapse:syn26320732
'''




# large_txt_to_mtx('./ica_bone_marrow_h5-filtered.txt','./data',gene_is_index=True)
# gene_ensembl = pd.read_csv('./data/genes.tsv',sep='\t',header=None).squeeze().tolist()
# gene_symbol = pd.Series(GeneConvert.ensemblgene_to_symbol(gene_ensembl,'human'))
# gene_symbol.to_csv('./data/genes_new.tsv',sep='\t',header=None,index=None)

# analyze big
big = pd.read_csv('big.txt',sep='\t',index_col=0)
adata = mtx_to_adata('./data')
adata = adata[big.index.tolist(),:]
adata = add_annotations(adata,'big.txt',cols=['Hay-et-al','ICGS2','Seurat3-Integration','Monocle3_cluster.txt'])

big_umap = pd.read_csv('ICGS2-Fig3-SI-UMAP-coordinates.txt',sep='\t',index_col=0,header=None)
big_umap_index = big_umap.index
big_umap_index = [item.split('.')[0] for item in big_umap_index]
big_umap.index = big_umap_index
big_umap.columns = ['umap_x','umap_y']
big_umap.to_csv('big_umap_trim.txt',sep='\t')

adata = add_umap(adata,'big_umap_trim.txt',mode='pandas',cols=['umap_x','umap_y'])
adata.obs.to_csv('check_obs.txt',sep='\t')
adata.write('./big_adata.h5ad')



# run sctriangulate using new version
adata = sc.read('big_adata.h5ad')
adata = just_log_norm(adata)
sctri = ScTriangulate(dir='./output_big_tfidf5_new',adata=adata,query=['Hay-et-al','ICGS2','Seurat3-Integration','Monocle3_cluster.txt'])
sctri.lazy_run(compute_metrics_parallel=False,scale_sccaf=False,win_fraction_cutoff=0.3)

sctri = ScTriangulate.deserialize('output_big_tfidf5_new/after_pruned_assess.p')
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')

'''DC'''
sctri.plot_heterogeneity('Hay-et-al','Dendritic Cell','umap',subset=['ICGS2@19','Seurat3-Integration@13','Seurat3-Integration@53',
                         'Seurat3-Integration@18','Seurat3-Integration@51'],merge=[('ICGS2@19','Seurat3-Integration@13')],format='png')
sctri.plot_heterogeneity('Hay-et-al','Dendritic Cell','heatmap',subset=['ICGS2@19','Seurat3-Integration@13','Seurat3-Integration@53',
                         'Seurat3-Integration@18','Seurat3-Integration@51'],merge=[('ICGS2@19','Seurat3-Integration@13')])
marker_gene_dict = {
    'Seurat3-Integration@18':['JCHAIN','MZB1','GZMB','IRF7','LILRA4'],
    'Seurat3-Integration@53':['PPP1R14A','AXL'],
    'ICGS2@19+Seurat3-Integration@13':['CLEC10A','CSTA','CEBPD'],
    'Seurat3-Integration@51':['C1orf54']
}
sctri.plot_heterogeneity('Hay-et-al','Dendritic Cell','heatmap+umap',subset=['ICGS2@19','Seurat3-Integration@13','Seurat3-Integration@53',
                         'Seurat3-Integration@18','Seurat3-Integration@51'],merge=[('ICGS2@19','Seurat3-Integration@13')],
                         marker_gene_dict=marker_gene_dict)

sctri.plot_heterogeneity('Hay-et-al','Dendritic Cell','multi_gene',subset=['ICGS2@19','Seurat3-Integration@13','Seurat3-Integration@53',
                         'Seurat3-Integration@18','Seurat3-Integration@51'],merge=[('ICGS2@19','Seurat3-Integration@13')],
                         multi_gene=['MZB1','AXL','CEBPD','C1orf54'])

'''Follicular B cell'''
sctri.plot_heterogeneity('Hay-et-al','Follicular B cell','umap',subset=['Monocle3_cluster.txt@11','Seurat3-Integration@50','Seurat3-Integration@18'],format='png')
sctri.plot_heterogeneity('Hay-et-al','Follicular B cell','heatmap',subset=['Monocle3_cluster.txt@11','Seurat3-Integration@50','Seurat3-Integration@18'])
marker_gene_dict = {
  'Seurat3-Integration@18':['CST3','C12orf75','SEC61B','ITM2C','TYROBP','APP'],
  'Seurat3-Integration@50':['IL32','CD3D','CD3E','CCL5','TRAC','TRBC1'],
  'Monocle3_cluster.txt@11':['MS4A1','CD79A','CD79B','BANK1','LINC01781','BLK','IGHA1']
}
sctri.plot_heterogeneity('Hay-et-al','Follicular B cell','heatmap+umap',subset=['Monocle3_cluster.txt@11','Seurat3-Integration@50','Seurat3-Integration@18'],
                            marker_gene_dict=marker_gene_dict)
sctri.plot_heterogeneity('Hay-et-al','Follicular_B_cell','multi_gene',subset=['Monocle3_cluster.txt@11','Seurat3-Integration@50','Seurat3-Integration@18'],multi_gene=['MS4A1','ITM2C','CD3D'])
sys.exit('stop')

'''Naive T-cell'''
sctri.plot_heterogeneity('Hay-et-al','Naive T-cell','umap',subset=['Hay-et-al@CD8 T-cell','Seurat3-Integration@4','Hay-et-al@Naive T-cell'],format='png')
sctri.plot_heterogeneity('Hay-et-al','Naive T-cell','heatmap',subset=['Hay-et-al@CD8 T-cell','Seurat3-Integration@4','Hay-et-al@Naive T-cell'])

marker_gene_dict = {
    'Seurat3-Integration@4':['CD8B','CD8A','LINC02446'],
    'Hay-et-al@CD8 T-cell':['CCL5','NKG7','CST7','GZMA','GZMK','DUSP2'],
    'Hay-et-al@Naive T-cell':['GPR183','SOCS3','ITGB1','CD4']
}
sctri.plot_heterogeneity('Hay-et-al','Naive T-cell','heatmap+umap',subset=['Hay-et-al@CD8 T-cell','Seurat3-Integration@4','Hay-et-al@Naive T-cell'],marker_gene_dict=marker_gene_dict)
sctri.plot_heterogeneity('Hay-et-al','Naive T-cell','multi_gene',subset=['Hay-et-al@CD8 T-cell','Seurat3-Integration@4','Hay-et-al@Naive T-cell'],multi_gene=['CD8A','NKG7','ITGB1'])



















# adata = sc.read('big_adata.h5ad')
# adata = just_log_norm(adata)
# sctri = ScTriangulate(dir='./output_big_tfidf5',adata=adata,query=['Hay-et-al','ICGS2','Seurat3-Integration','Monocle3_cluster.txt'],
#                         reference='Hay-et-al',species='human',criterion=2,verbose=2)
# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})
# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=False)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley()
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri = ScTriangulate.deserialize('output_big_tfidf5/after_shapley.p')
# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_prune_rank.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_big_tfidf5/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()

# sctri.add_to_invalid_by_win_fraction(percent=0.3)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Hay-et-al')
# sctri.plot_umap('confidence','continuous',umap_cmap='viridis')
# for col in ['final_annotation','raw','pruned']:
#     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='Hay-et-al')
# sctri.viewer_heterogeneity_html(key='Hay-et-al')
# sys.exit('stop')


















# adata = sc.read('small_adata.h5ad')
# adata = just_log_norm(adata)
# sctri = ScTriangulate(dir='./output_small_tfidf5_no_raceid',adata=adata,query=['Hay-et-all','ICGS2','Seurat3','Monocle2','SC3','Monocle3','Scanpy'],
#                         reference='Hay-et-all',species='human',criterion=2,verbose=2)

# sctri.add_new_metrics({'tfidf5':tf_idf5_for_cluster})

# sctri.doublet_predict()  
# sctri.compute_metrics(parallel=True)
# sctri.serialize(name='after_metrics.p')
# sctri.obs_to_df(name='sctri_obs_after_metrics.txt')

# sctri.compute_shapley(parallel=True)
# sctri.serialize(name='after_shapley.p')
# sctri.obs_to_df(name='sctri_obs_after_shapley.txt')

# sctri.pruning(method='rank',discard=None)
# sctri.serialize(name='after_prune_rank.p')
# sctri.uns['raw_cluster_goodness'].to_csv('output_small_tfidf5_no_raceid/raw_cluster_goodness.txt',sep='\t')

# sctri.viewer_cluster_feature_figure(parallel=False)
# sctri.viewer_cluster_feature_html()


# sctri = ScTriangulate.deserialize('output_small_tfidf5_no_raceid/after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.3)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Hay-et-all')
# for col in ['final_annotation','raw','pruned']:
#     sctri.plot_umap(col,'category',True,'pdf')
# sctri.viewer_heterogeneity_figure(key='Hay-et-all')
# sctri.viewer_heterogeneity_html(key='Hay-et-all')


















