#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os
import sys
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import re
from sctriangulate import ScTriangulate
from sctriangulate.colors import bg_greyed_cmap
from sctriangulate.preprocessing import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

'''
all the input and intermediate files are on https://www.synapse.org/#!Synapse:syn26320419
'''

adata = sc.read_10x_h5('pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5',gex_only=False)  # 11909
adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
adata_atac = adata[:,adata.var['feature_types']=='Peaks']
adata_atac = reformat_peak(adata_atac)
find_genes(adata_atac,'gencode.v38.annotation.gtf') 
adata_atac.var_names = [str(v) + '_' + adata_atac.var['gene_annotation'][i] for i,v in enumerate(adata_atac.var_names)]

# QC, RNA
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
    sc.pl.violin(adata_rna,key,jitter=0.4)
    plt.savefig('qc_rna_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_rna,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
plt.savefig('qc_rna_scatter.pdf',bbox_inches='tight')
plt.close()

# QC, ATAC
sc.pp.calculate_qc_metrics(adata_atac, var_type='peaks',percent_top=None, log1p=False, inplace=True)
for key in ['n_peaks_by_counts','total_counts']:
    sc.pl.violin(adata_atac,key,jitter=0.4)
    plt.savefig('qc_atac_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_atac,x='n_peaks_by_counts',y='total_counts')
plt.savefig('qc_peak_scatter.pdf',bbox_inches='tight')
plt.close()

# decision
for RNA, min_gene = 300, min_counts = 500, mito < 20%
for ATAC, min_peak = 1000
adata_rna = adata_rna[adata_rna.obs['n_genes_by_counts']>300,:]
adata_rna = adata_rna[adata_rna.obs['total_counts']>500,:]
adata_rna = adata_rna[adata_rna.obs['pct_counts_mt']<20,:]
kept_rna = adata_rna.obs_names
adata_atac = adata_atac[adata_atac.obs['n_peaks_by_counts']>1000,:]
kept_atac = adata_atac.obs_names
common = list(set(kept_rna).intersection(set(kept_atac)))

adata_rna = adata_rna[adata_rna.obs_names.isin(common),:]   # 10991 × 36601
adata_atac = adata_atac[adata_atac.obs_names.isin(common),:]  # 10991 × 108344
adata_rna.write('post_rna_qc.h5ad')  # this file will serve as input to Azimuth web app

# get the clusters and combine
doublet_map = doublet_predict(adata_rna)
adding_azimuth(adata_rna,'azimuth_pred.tsv')
adata_rna = scanpy_recipe(adata_rna,False,resolutions=[0.5,1,2,3,4,5,6],modality='rna',pca_n_comps=50,n_top_genes=3000)
adata_atac = scanpy_recipe(adata_atac,False,resolutions=[0.5,1,2,3,4,5,6],modality='atac',pca_n_comps=100,n_top_genes=100000)

adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_atac = sc.read('adata_after_scanpy_recipe_atac_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_combine = concat_rna_and_other(adata_rna,adata_atac,umap='other',name='atac',prefix='')
adata_combine = make_sure_adata_writable(adata_combine,delete=True)
adata_combine.write('combined_rna_atac_au.h5ad')

# plot
all_r = ['sctri_{}_leiden_{}'.format(m,i) for i in [0.5,1,2,3,4,5,6] for m in ['rna','atac']]
cols = all_r + ['azimuth','doublet_scores']
umap_dual_view_save(adata_combine,cols)


# run scTriangulate
adata_combine = sc.read('combined_rna_atac_ru.h5ad')
sctri = ScTriangulate(dir='output_two',adata=adata_combine,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3',
                      'sctri_atac_leiden_1','sctri_atac_leiden_2','sctri_atac_leiden_3'],predict_doublet='precomputed')
sctri.lazy_run(win_fraction_cutoff=0.35,scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_atac_leiden_1'])


## insights
sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')
sctri.modality_contributions()
for col in ['rna_contribution','atac_contribution']:
   sctri.plot_umap(col,'continuous',umap_cmap='viridis')
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')


# CD4 populations
for cluster in ['sctri_rna_leiden_3@5','sctri_rna_leiden_2@14','sctri_rna_leiden_1@13','sctri_rna_leiden_3@24','sctri_rna_leiden_3@23','sctri_rna_leiden_3@1','sctri_rna_leiden_2@21']:
   sctri.plot_multi_modal_feature_rank(cluster=cluster)

# let's focus on the CD4 population and make some figures
sctri.adata.obs['dummy_key'] = np.full(shape=sctri.adata.obs.shape[0],fill_value='dummy_cluster')
subset = ['sctri_rna_leiden_2@21','sctri_rna_leiden_3@1','sctri_rna_leiden_3@23','sctri_rna_leiden_3@5','sctri_rna_leiden_3@24','sctri_rna_leiden_2@14','sctri_rna_leiden_1@13']
sctri.plot_heterogeneity('dummy_key','dummy_cluster','umap','pruned',subset=subset)
for gene in ['PTPN13','SLC4A10','IL2RA','ITGB1','FHIT','TSHZ2','chr20_52971902_52974673_TSHZ2','chr6_22045235_22045691_CASC15','chr6_22146145_22147503_NBAT1;CASC15','SOX4','chr6_22148300_22148794_NBAT1;CASC15','chr6_21876998_21877505_CASC15','CASC15']:
   sctri.plot_heterogeneity('dummy_key','dummy_cluster','single_gene','pruned',subset=subset,single_gene=gene,cmap='viridis')
sys.exit('stop')
subset = ['CD4 Naive','CD4 TCM','Treg','CD4 TEM','MAIT']
sctri.plot_heterogeneity('dummy_key','dummy_cluster','umap','azimuth',subset=subset)

sctri.plot_umap('CASC15','continuous',umap_cmap='viridis')
sctri.plot_umap('MAP3K4','continuous',umap_cmap='viridis')
genes = ['CXCR6','GZMB','PRDM1','CXCR5','IL21','BCL6','CCR7','IL7R','S1PR1']
sc.pl.umap(sctri.adata,color=genes,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('cd4_markers.pdf',bbox_inches='tight')
plt.close()






























































