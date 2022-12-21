import os
import sys
from sctriangulate import *
from sctriangulate.preprocessing import *

# harmony_sct
adata = sc.read('harmony_sct.1.h5ad')

'''
adata
AnnData object with n_obs × n_vars = 363356 × 32284
    obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'DonorID', 'DataID'
    var: '_index', 'features'
    obsm: 'X_pca', 'X_umap'


adata.raw.X
<363356x32284 sparse matrix of type '<class 'numpy.float64'>'
	with 822807723 stored elements in Compressed Sparse Row format>

adata.raw.X[:20,:20].toarray()  # confirm they are raw count
adata.X[:20,:20].toarray()      # confirm they are integrated data
adata.X.min() return 0, means it is non-negative, no need to supply another copy in layer



My comments:

let me know if I understand correctly, you first run sctransform on all the data, and run harmony on sctransform-derived pca space to get the harmonized pca space.
Then you would like to derive clusters in different resolutions and assess the stability respectively. 

If that's the case, just to confirm, is the "X_pca" shown in obsm the post-harmony PCA space? If so, then very easy, you first need to get the clusters using following code

sc.pp.neighbors(adata,use_rep='X_pca')
sc.tl.leiden(adata,resolution=1,key_added='sctri_{}_leiden_{}'.format('rna',1))  # you can get multiple resolution by writing for loop

Then run scTriangulate

del adata.raw
sctri = ScTriangulate(dir='./output_harmony_sct',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
sctri.compute_metrics(parallel=False)
'''


# mnn_sct
adata = sc.read('mnn_sct.1.h5ad')

'''
>>> adata
AnnData object with n_obs × n_vars = 363356 × 32284
    obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'DonorID', 'DataID'
    var: '_index', 'features'
    obsm: 'X_pca'


adata.raw.X
<363356x32284 sparse matrix of type '<class 'numpy.float64'>'
	with 822807723 stored elements in Compressed Sparse Row format>

adata.raw.X[:20,:20].toarray()  # confirm they are raw count
adata.X[:20,:20].toarray()      # confirm they are integrated data
adata.X.min() return 0, means it is non-negative, no need to supply another copy in layer


My comments:

If you first run sct, and then mnn, just confirm how the 'X_pca' is derived, if it is pca space derived after MNN, then you just follow exact same comment I wrote for harmony. 
'''


# rpca_sct
'''
I think it's the same instruction as the harmony one

'''






