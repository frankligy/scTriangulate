import pytest
import os

os.chdir(os.path.dirname(__file__))

import scanpy as sc
from sctriangulate import *
from sctriangulate.colors import *
from sctriangulate.preprocessing import *
from sctriangulate.spatial import *

'''
some useful pytest points, from https://www.tutorialspoint.com/pytest/index.htm

pytest mini_test.py -k import -m plot -v
@pytest.mark.plot
@pytest.fixture
@pytest.mark.parametrize
@pytest.mark.xfail
@pytest.mark.skip

'''

@pytest.mark.main
def test_global_setting():
    sctriangulate_setting(backend='Agg')

@pytest.mark.main
def test_main_function():
    adata = sc.read('input.h5ad')
    sctri = ScTriangulate(dir='output',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
    sctri.lazy_run()

@pytest.mark.output
def test_display_hierarchy():
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.add_to_invalid_by_win_fraction(percent=0.25)
    sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=sctri.reference)
    sctri.display_hierarchy(ref_col='sctri_rna_leiden_1',query_col='raw')

@pytest.mark.parametrize('style',['umap'])
@pytest.mark.plot
def test_plot_heterogeneity(style):
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='3',col='pruned',style=style)

@pytest.mark.parametrize('cluster, feature',[('3','enrichment'),('3','marker_genes'),('3','exclusive_genes'),('3','location')])
@pytest.mark.plot
def test_plot_cluster_feature(cluster,feature):
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_cluster_feature(key='sctri_rna_leiden_1',cluster=cluster,feature=feature)

@pytest.mark.plot
def test_plot_winners_statistics():
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_winners_statistics(col='raw',fontsize=6)

@pytest.mark.plot
def test_plot_clusterability():
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_clusterability(key='sctri_rna_leiden_1',col='raw',fontsize=8)

@pytest.mark.plot
def test_plot_umap():
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_umap(col='pruned',kind='category')
    sctri.plot_umap(col='confidence',kind='continuous',umap_cmap='viridis')

@pytest.mark.plot
def test_plot_long_heatmap(): 
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')  
    sctri.plot_long_heatmap(key='sctri_rna_leiden_1',n_features=20,figsize=(20,20))

@pytest.mark.plot
def test_plot_confusion():
    sctriangulate_setting(backend='Agg')
    sctri = ScTriangulate.deserialize('output/after_rank_pruning.p')
    sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1',cmap=retrieve_pretty_cmap('shap'))

























