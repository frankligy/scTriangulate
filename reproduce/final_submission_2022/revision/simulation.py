#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
import os,sys
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# poc
adata = sc.read('sim_poc.h5ad')
'''
AnnData object with n_obs × n_vars = 3000 × 10000
    obs: 'Cell', 'Batch', 'Group', 'ExpLibSize', 'sizeFactor'
    var: 'Gene', 'BaseGeneMean', 'OutlierFactor', 'GeneMean', 'DEFacGroup1', 'DEFacGroup2', 'DEFacGroup3', 'DEFacGroup4', 'DEFacGroup5'
    uns: 'X_name'
    obsm: 'PCA'
    layers: 'BCV', 'BaseCellMeans', 'CellMeans', 'TrueCounts', 'counts', 'logcounts'
'''


adata.X = adata.layers['counts']
adata.obsm['X_umap'] = adata.obsm['PCA']

adata.obs['broad_anno'] = adata.obs['Group'].map({'Group5':'c1','Group4':'c2','Group2':'c3','Group1':'c3','Group3':'c3'}).values
adata.obs['middle_anno'] = adata.obs['Group'].map({'Group5':'c1','Group4':'c2','Group2':'c3','Group1':'c4','Group3':'c3'}).values
adata.obs['fine_anno'] = adata.obs['Group'].map({'Group5':'c1','Group4':'c2','Group2':'c3','Group1':'c4','Group3':'c5'}).values

col = []
for item in adata.obs['Group']:
    if item != 'Group2' and item != 'Group3':
        col.append(item)
    else:
        suffix = str(np.random.choice([1,2],size=1)[0])
        col.append(item + '+' + suffix)
adata.obs['overcluster'] = col

umap_dual_view_save(adata,cols=['broad_anno','middle_anno','fine_anno','overcluster','Group'])

sctri = ScTriangulate(dir='output_sim_poc_one_overcluster',adata=adata,query=['broad_anno','middle_anno','fine_anno','overcluster'],add_metrics={})
sctri.lazy_run(scale_sccaf=False,added_metrics_kwargs={},assess_pruned=True,viewer_cluster=False,viewer_heterogeneity=False)
sys.exit('stop')


# # discovery mode and conserved mode
# adata = sc.read('sim_dc.h5ad')
# adata.X = adata.layers['counts']
# adata.obsm['X_umap'] = adata.obsm['PCA']
# adata.obs['broad_anno'] = adata.obs['Group'].map({'Group1':'c1','Group2':'c2','Group3':'c3','Group4':'c3','Group5':'c3','Group6':'c3'}).values
# adata.obs['middle_anno'] = adata.obs['Group'].map({'Group1':'c1','Group2':'c2','Group3':'c3','Group4':'c3','Group5':'c4','Group6':'c4'}).values
# adata.obs['fine_anno'] = adata.obs['Group'].map({'Group1':'c1','Group2':'c2','Group3':'c3','Group4':'c4','Group5':'c5','Group6':'c6'}).values
# umap_dual_view_save(adata,cols=['broad_anno','middle_anno','fine_anno','Group'])
# sctri = ScTriangulate(dir='output_sim_dc_one',adata=adata,query=['broad_anno','middle_anno','fine_anno'],add_metrics={})
# sctri.lazy_run(scale_sccaf=False,added_metrics_kwargs={},assess_pruned=True,viewer_cluster=False,viewer_heterogeneity=False)
# sctri = ScTriangulate(dir='output_sim_dc_two',adata=adata,query=['broad_anno','middle_anno','fine_anno'])
# sctri.lazy_run(scale_sccaf=False,assess_pruned=True,viewer_cluster=False,viewer_heterogeneity=False)