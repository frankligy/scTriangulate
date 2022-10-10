import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import scvi
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import seaborn as sns
import scanpy as sc
import anndata as ad
import os,sys

os.chdir('.')



# Load the data for opentarget
adata_count = sc.read('cell2loc_input_adata_sc.h5ad')
adata_st1 = sc.read('cell2loc_input_adata_spatial.h5ad')

# remove mitochondrial genes
adata_count.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_count.var_names]
adata_count.obsm['MT'] = adata_count[:, adata_count.var['MT_gene'].values].X.toarray()
adata_count = adata_count[:, ~adata_count.var['MT_gene'].values]
adata_st1.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_st1.var_names]
adata_st1.obsm['MT'] = adata_st1[:, adata_st1.var['MT_gene'].values].X.toarray()
adata_st1 = adata_st1[:, ~adata_st1.var['MT_gene'].values]

# get signature reference
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_count, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_count = adata_count[:,selected]
from cell2location.models import RegressionModel
RegressionModel.setup_anndata(adata=adata_count,
                        batch_key='sample',
                        labels_key='final_annotation')
mod = RegressionModel(adata_count)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)
mod.plot_history(20)
adata_count = mod.export_posterior(adata_count, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})
mod.save('regression_model', overwrite=True)
adata_count.write('regression_ref.h5ad')
mod.plot_QC()

# cell2location
mod = RegressionModel.load('regression_model', adata_count)
adata_count = sc.read_h5ad('regression_ref.h5ad')
inf_aver = adata_count.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_count.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_count.uns['mod']['factor_names']
intersect = np.intersect1d(adata_st1.var_names, inf_aver.index)
adata_st1 = adata_st1[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
from cell2location.models import Cell2location
adata_st1.obs['sample'] = np.full(adata_st1.shape[0],fill_value='sample1')
Cell2location.setup_anndata(adata=adata_st1, batch_key="sample")
mod = Cell2location(adata_st1, cell_state_df=inf_aver,N_cells_per_location=30,detection_alpha=200)
mod.train(max_epochs=30000,batch_size=None,train_size=1,use_gpu=True)
mod.plot_history(1000)
adata_st1 = mod.export_posterior(adata_st1, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True})
mod.save('cell2location_model', overwrite=True)
adata_st1.write('cell2location_st1.h5ad')

# extract abundance
mod = Cell2location.load('cell2location_model', adata_st1)
adata_st1 = sc.read_h5ad('cell2location_st1.h5ad')
adata_st1.obsm['q05_cell_abundance_w_sf'].to_csv('cell2location_prop.txt',sep='\t')




