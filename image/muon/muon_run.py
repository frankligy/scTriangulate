#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/Hs_BM_400k_RNA/muon_env/bin/python3.8

import muon as mu
from muon import MuData
import numpy as np
import scanpy as sc
import pandas as pd


import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# adata_rna = sc.read('to_muon_rna.h5ad')
# adata_adt = sc.read('to_muon_adt.h5ad')

# # re-derive neighbors
# sc.pp.neighbors(adata_rna,use_rep='X_pca_harmony')
# sc.pp.neighbors(adata_adt,use_rep='X_pca')


# # create mdata and do wnn
# mdata = MuData({'rna':adata_rna,'adt':adata_adt})
# mu.pp.neighbors(mdata)
# mdata.write('mudata.h5mu')

# mdata = mu.read('mudata.h5mu')
# mu.tl.umap(mdata)
# for r in [0.5,1,2,3,4,5,6]:
#     sc.tl.leiden(mdata, resolution=r,key_added='leiden_wnn_{}'.format(r))
# mdata.write('mudata_umap_leiden.h5mu')

mdata = mu.read('mudata_umap_leiden.h5mu')
# for key in ['leiden_wnn_{}'.format(r) for r in [0.5,1,2,3,4,5,6]]:
#     mu.pl.umap(mdata, color=[key],legend_loc='on data',legend_fontsize='xx-small')
#     plt.savefig('muon_wnn_umap_{}.pdf'.format(key.replace(':','_')),bbox_inches='tight')
#     plt.close()

pd.DataFrame(data=mdata.obsm['X_umap'],index=mdata.obs_names,columns=['umap_x','umap_y']).to_csv('muon_wnn_umap.txt',sep='\t')
mdata.obs.to_csv('muon_wnn_metadata.txt',sep='\t')
sys.exit('stop')

mu.pl.umap(mdata, color=['rna:mod_weight', 'adt:mod_weight'])
plt.savefig('muon_wnn_umap_weights.pdf',bbox_inches='tight')
plt.close()
for key in ['rna:ori_label','rna:sample','rna:donor']:
    mu.pl.umap(mdata, color=[key],legend_loc='on data',legend_fontsize='xx-small')
    plt.savefig('muon_wnn_umap_{}.pdf'.format(key.replace(':','_')),bbox_inches='tight')
    plt.close()





