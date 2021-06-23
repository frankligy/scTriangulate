
'''This script is to test sctriangulate program,
   Using pbmc3k dataset'''

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
os.chdir('/Users/ligk2e/Desktop/scTriangulate')
sys.path.append('.')
from sctriangulate import *

# load the data
adata = sc.read('pbmc3k_azimuth_umap.h5ad')

# instantiation
sctri = ScTriangulate(dir='./output',adata=adata,query=['leiden1','leiden2','leiden3'])

# main program
sctri.compute_metrics(parallel=True,scale_sccaf=False)
sctri.compute_shapley(parallel=True)
sctri.pruning(method='rank',scale_sccaf=False)

# clean step
sctri.add_to_invalid_by_win_fraction(percent=0.25)
sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='leiden1')

# IO
sctri.serialize(name='after_prune_rank.p')
sctri = ScTriangulate.deserialize('output/after_prune_rank.p')


# generate viewer
sctri.viewer_cluster_feature_figure()
sctri.viewer_cluster_feature_html()
sctri.add_to_invalid_by_win_fraction(percent=0.25)
sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='leiden3')
sctri.viewer_heterogeneity_figure(keys=['leiden3'])
sctri.viewer_heterogeneity_html(key='leiden3')

'''test plot_heterogeneity function'''
# umap
sctri.plot_heterogeneity('leiden1','0','umap',subset=['leiden1@0','leiden3@10'])

# heatmap
sctri.plot_heterogeneity('leiden1','0','heatmap',subset=['leiden1@0','leiden3@10'])

# heatmap_custom_gene
marker_gene_dict = {
    'leiden1@0':['MAPK14','RASSF1','PARP10'],
    'leiden3@10':['ANXA1','CD52','CMTM7']
}
sctri.plot_heterogeneity('leiden1','0','heatmap_custom_gene',subset=['leiden1@0','leiden3@10'],marker_gene_dict=marker_gene_dict)

# violin
sctri.plot_heterogeneity('leiden1','0','violin',subset=['leiden1@0','leiden3@10'],genes=['MAPK14','ANXA1'])

# sankey
sctri.plot_heterogeneity('leiden1','0','sankey')

# cellxgene
sctri.plot_heterogeneity('leiden1','0','cellxgene')

# build
sctri.plot_heterogeneity('leiden1','0','build')

'''test other plotting function'''
sctri.plot_circular_barplot('leiden1','pruned')
sctri.plot_confusion('confusion_reassign','leiden1',annot=True)
sctri.plot_cluster_feature('leiden1','0','enrichment')
sctri.plot_cluster_feature('leiden1','0','marker_genes')
sctri.plot_cluster_feature('leiden1','0','exclusive_genes')
sctri.plot_cluster_feature('leiden1','0','location')
sctri.plot_umap('confidence','continuous',umap_dot_size=10)
df = sctri.plot_winners_statistics('raw')
bucket = sctri.plot_clusterability('pruned')

# output useful intermediate result
df = sctri.get_metrics_and_shapley('TTTCTACTGAGGCA-1')
sctri.obs_to_df()
sctri.var_to_df()
sctri.gene_to_df('exclusive_genes','leiden1')
sctri.confusion_to_df('confusion_sccaf','leiden1')
sctri.display_hierarchy('pruned',True)























