
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
os.chdir('/Users/ligk2e/Desktop/scTriangulate')
sys.path.append('.')

from sctriangulate import *


adata = sc.read('pbmc3k_azimuth_umap.h5ad')
sctri = ScTriangulate(dir='./output',adata=adata,query=['leiden1','leiden2','leiden3'])

sctri.compute_metrics(parallel=True,scale_sccaf=False)
# sctri.penalize_artifact(mode='void',stamps=['leiden1@2'])

sctri.compute_shapley(parallel=True)

sctri.pruning(method='rank',scale_sccaf=False)

sctri.viewer_cluster_feature_figure()
sctri.viewer_cluster_feature_html()

sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='leiden3')
sctri.viewer_heterogeneity_figure(keys=['leiden3'])
sctri.viewer_heterogeneity_html(key='leiden3')

sctri.get_cluster()





