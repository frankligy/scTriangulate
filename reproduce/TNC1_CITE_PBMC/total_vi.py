#! /data/salomonis2/LabFiles/Frank-Li/scTriangulate/TNC1_qc/scvi_env/bin/python3.7

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scvi
import scanpy as sc
import pickle
import os,sys



# # tutorial running
# adata = scvi.data.pbmcs_10x_cite_seq(run_setup_anndata=False)
# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# adata.raw = adata
# sc.pp.highly_variable_genes(adata,n_top_genes=4000,flavor="seurat_v3",batch_key="batch",subset=True,layer="counts")
# scvi.model.TOTALVI.setup_anndata(adata,protein_expression_obsm_key="protein_expression",layer="counts",batch_key="batch")
# vae = scvi.model.TOTALVI(adata, latent_distribution="normal")
# print(vae.train())

# plt.plot(vae.history["elbo_train"], label="train")
# plt.plot(vae.history["elbo_validation"], label="validation")
# plt.title("Negative ELBO over training epochs")
# plt.ylim(1200, 1400)
# plt.legend()
# plt.savefig('check_train.pdf',bbox_inches='tight')
# plt.close()

# adata.obsm["X_totalVI"] = vae.get_latent_representation()
# sc.pp.neighbors(adata, use_rep="X_totalVI")
# sc.tl.umap(adata, min_dist=0.4)
# sc.tl.leiden(adata, key_added="leiden_totalVI")
# sc.pl.umap(adata,color=["leiden_totalVI", "batch"],frameon=False,ncols=1)
# plt.savefig('check_umap.pdf',bbox_inches='tight')
# plt.close()

# own running
adata = sc.read_10x_h5('28WM_ND19-341__TNC-RNA-ADT.h5',gex_only=False)
adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
adata_adt = adata[:,adata.var['feature_types']=='Antibody Capture']  # 8491
adata_rna.var_names_make_unique()
adata_adt.var_names_make_unique()
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata_rna, min_genes=300)
sc.pp.filter_cells(adata_rna, min_counts=500)
adata_rna = adata_rna[adata_rna.obs.pct_counts_mt < 20, :]
adata_adt = adata_adt[adata_rna.obs_names,:]   # 6406

adata_scvi = adata_rna.copy()
adata_scvi.obsm['protein_expression'] = adata_adt.to_df()
adata_scvi.layers['counts'] = adata_scvi.X.copy()
sc.pp.normalize_total(adata_scvi, target_sum=1e4)
sc.pp.log1p(adata_scvi)
adata_scvi.raw = adata_scvi
sc.pp.highly_variable_genes(adata_scvi,n_top_genes=4000,flavor="seurat_v3",subset=True,layer="counts")
scvi.model.TOTALVI.setup_anndata(adata_scvi,protein_expression_obsm_key="protein_expression",layer="counts")
vae = scvi.model.TOTALVI(adata_scvi, latent_distribution="normal")
print(vae.train())

plt.plot(vae.history["elbo_train"], label="train")
plt.plot(vae.history["elbo_validation"], label="validation")
plt.title("Negative ELBO over training epochs")
#plt.ylim(1200, 1400)
plt.legend()
plt.savefig('check_own_train.pdf',bbox_inches='tight')
plt.close()

adata_scvi.obsm["X_totalVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata_scvi, use_rep="X_totalVI")
sc.tl.umap(adata_scvi, min_dist=0.4)
for r in [0.5,1,2,3,4,5,6]:
    sc.tl.leiden(adata_scvi,resolution=r,key_added="leiden_totalVI_{}".format(r))
    sc.pl.umap(adata_scvi,color="leiden_totalVI_{}".format(r),frameon=False,ncols=1)
    plt.savefig('check_own_umap_{}.pdf'.format(r),bbox_inches='tight')
    plt.close()
adata_scvi.obs.to_csv('totalVI_obs.txt',sep='\t')
