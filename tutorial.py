
import os
import sys
import scanpy as sc
os.chdir('/Users/ligk2e/Desktop/scTriangulate')
sys.path.append('.')
sys.path.append('./sctriangulate')

from sctriangulate import ScTriangulate

adata = sc.read('pbmc3k_azimuth_umap.h5ad')
sctri = ScTriangulate(dir='.',adata=adata,query=['leiden1','leiden2','leiden3'],reference='leiden1',
                      species='human',criterion=2)
sctri.doublet_predict()  # after this step. X is dense
# sctri.plot_umap('doublet_scores','continuous',True)
# sctri.plot_umap('leiden1','category',True)

sctri.compute_metrics()
# sctri.plot_confusion('confusion_sccaf','leiden3',True)
# sctri.plot_cluster_feature('leiden2','3','exclusive_gene',True)

sctri.compute_shapley()
# sctri.plot_umap('raw','category',False)

sctri.pruning()
# sctri.plot_umap('pruned','category',False)
# sctri.plot_umap('prefixed','category',False)
# sctri.plot_heterogeneity('1','umap',True)


sctri.building_viewer_fig()
sctri.building_viewer_html()

sctri.get_cluster()
# sctri.plot_umap('user_choice','continuous',False)




