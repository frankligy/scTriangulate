import scanpy as sc
from sctriangulate import ScTriangulate
from sctriangulate.colors import *
from sctriangulate.preprocessing import *

adata = sc.read('input.h5ad')
sctri = ScTriangulate(dir='output_test',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],predict_doublet=False)
sctri.lazy_run(assess_pruned=True)



