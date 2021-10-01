import scanpy as sc
from sctriangulate import *
from sctriangulate.colors import *
from sctriangulate.preprocessing import *
print('===================\nimport modules test:',u'\u2713','\n====================')

sctriangulate_setting(backend='Agg')

adata = sc.read('input.h5ad')
sctri = ScTriangulate(dir='output',adata=adata,add_metrics={},query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
print('====================\ninstantiation test:',u'\u2713','\n====================')

sctri.lazy_run(scale_sccaf=False,viewer_cluster=False,viewer_heterogeneity=False)
print('=====================\nlazy_run test:',u'\u2713','\n========================')

sctri.plot_winners_statistics(col='raw',fontsize=6)
print('======================\nplot_winners_statistics test:',u'\u2713','\n=======================')
sctri.plot_clusterability(key='sctri_rna_leiden_1',col='raw',fontsize=8)
print('======================\nplot_clusterability test:',u'\u2713','\n=========================')
sctri.display_hierarchy(ref_col='sctri_rna_leiden_1',query_col='raw')
print('======================\ndisplay_hierarchy test:',u'\u2713','\n=====================')
sctri.plot_umap(col='pruned',kind='category')
sctri.plot_umap(col='confidence',kind='continuous',umap_cmap='viridis')
print('======================\nplot_umap test:',u'\u2713','\n===================')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='enrichment')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='marker_genes')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='exclusive_genes')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='location')
print('======================\nplot_cluster_feature test:',u'\u2713','\n=====================')
sctri.plot_long_heatmap(key='pruned',n_features=20,figsize=(20,20))
print('=======================\nplot_long_heatmap test:',u'\u2713','\n========================')
sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1',cmap=retrieve_pretty_cmap('shap'))
print('=======================\nplot_confusion test:',u'\u2713','\n============================')
























