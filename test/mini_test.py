import scanpy as sc
from sctriangulate import ScTriangulate
from sctriangulate.colors import *
from sctriangulate.preprocessing import *
print('import modules test:',u'\u2713')

adata = sc.read('input.h5ad')
sctri = ScTriangulate(dir='output',adata=adata,add_metrics={},query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
print('instantiation test:',u'\u2713')

sctri.lazy_run(scale_sccaf=False)
print('lazy_run test:',u'\u2713')

sctri.plot_winners_statistics(col='raw',fontsize=6)
print('plot_winners_statistics test:',u'\u2713')
sctri.plot_clusterability(key='sctri_rna_leiden_1',col='raw',fontsize=8)
print('plot_clusterability test:',u'\u2713')
sctri.display_hierarchy(ref_col='sctri_rna_leiden_1',query_col='raw')
print('display_hierarchy test:',u'\u2713')
sctri.plot_umap(col='pruned',kind='category')
sctri.plot_umap(col='confidence',kind='continuous',umap_cmap='viridis')
print('plot_umap test:',u'\u2713')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='enrichment')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='marker_genes')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='exclusive_genes')
sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster='3',feature='location')
print('plot_cluster_feature test:',u'\u2713')
sctri.plot_long_heatmap(key='pruned',n_features=20,figsize=(20,20))
print('plot_long_heatmap test:',u'\u2713')
sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1',cmap=retrieve_pretty_cmap('shap'))
print('plot_confusion test:',u'\u2713')
























