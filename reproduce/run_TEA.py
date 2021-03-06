#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os,sys

from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

'''
the input files are in https://www.synapse.org/#!Synapse:syn26320350
'''

# readin the original data from GEO, and process a bit
# the data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4949911
adata = sc.read_10x_h5('GSM4949911_X061-AP0C1W1_leukopak_perm-cells_tea_fulldepth_cellranger-arc_filtered_feature_bc_matrix.h5',gex_only=False)
adt = pd.read_csv('GSM4949911_X061-AP0C1W1_leukopak_perm-cells_tea_fulldepth_adt_counts.csv',sep=',',index_col=0)
# although adt count matrix processed by BarCounter, the index is unique, so let's append by "_1"
adt.index = [item + '-1' for item in adt.index]
adt_valid = adt.loc[adata.obs_names.tolist(),:]
adt_valid.drop(columns=['total'],inplace=True)
adt_valid.to_csv('ADT_count_matrix_processed.txt',sep='\t')

# get tri-modality adata
adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']   # 8213 × 36601
adata_atac = adata[:,adata.var['feature_types']=='Peaks']           # 8213 × 66828
adata_adt = small_txt_to_adata('ADT_count_matrix_processed.txt',False)  # 8213 × 46

# QC RNA
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
    sc.pl.violin(adata_rna,key,jitter=0.4)
    plt.savefig('qc_rna_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_rna,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
plt.savefig('qc_rna_scatter.pdf',bbox_inches='tight')
plt.close()

# QC ADT
sc.pp.calculate_qc_metrics(adata_adt, var_type='adts',percent_top=None, log1p=False, inplace=True)

for key in ['n_adts_by_counts','total_counts']:
    sc.pl.violin(adata_adt,key,jitter=0.4)
    plt.savefig('qc_adt_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_adt,x='n_adts_by_counts',y='total_counts')
plt.savefig('qc_adt_scatter.pdf',bbox_inches='tight')
plt.close()

# QC ATAC
sc.pp.calculate_qc_metrics(adata_atac, var_type='peaks',percent_top=None, log1p=False, inplace=True)
for key in ['n_peaks_by_counts','total_counts']:
    sc.pl.violin(adata_atac,key,jitter=0.4)
    plt.savefig('qc_atac_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_atac,x='n_peaks_by_counts',y='total_counts')
plt.savefig('qc_peak_scatter.pdf',bbox_inches='tight')
plt.close()


# decision, min_gene = 300, min_count = 500, mt < 20%, min_peak = 1000
adata_rna = adata_rna[adata_rna.obs['n_genes_by_counts']>300,:]
adata_rna = adata_rna[adata_rna.obs['total_counts']>500,:]
adata_rna = adata_rna[adata_rna.obs['pct_counts_mt']<20,:]
kept_rna = adata_rna.obs_names   # down to 7230
adata_atac = adata_atac[adata_atac.obs['n_peaks_by_counts']>1000,:]
kept_atac = adata_atac.obs_names   # down to 7639
common = list(set(kept_rna).intersection(set(kept_atac)))   # 6819

adata_rna = adata_rna[adata_rna.obs_names.isin(common),:]    # 6819 × 36601
adata_adt = adata_adt[adata_adt.obs_names.isin(common),:]    # 6819 × 66828
adata_atac = adata_atac[adata_atac.obs_names.isin(common),:]   # 6819 × 46
adata_rna.write('post_rna_qc.h5ad')  # this post-rna h5ad file will be subjected to Azimuth to get the mapping 

'''
The Azimuth prediction result named 'azimuth_pred.tsv' will be on https://www.synapse.org/#!Synapse:syn26320350
'''


# analyze each modality one by one
'''RNA'''
adding_azimuth(adata_rna,'azimuth_pred.tsv')
doublet_map = doublet_predict(adata_rna)
adata_rna = scanpy_recipe(adata_rna,False,resolutions=[0.5,1,2,3,4,5,6],modality='rna',pca_n_comps=50,n_top_genes=3000)

'''ADT'''
adding_azimuth(adata_adt,'azimuth_pred.tsv')
adata_adt.obs['doublet_scores'] = adata_adt.obs_names.map(doublet_map).values
adata_adt = scanpy_recipe(adata_adt,False,resolutions=[0.5,1,2,3,4,5,6],modality='adt',pca_n_comps=15)

'''ATAC'''
adding_azimuth(adata_atac,'azimuth_pred.tsv')
adata_atac.obs['doublet_scores'] = adata_atac.obs_names.map(doublet_map).values
adata_atac = reformat_peak(adata_atac,True)
find_genes(adata_atac,'gencode.v38.annotation.gtf')
adata_atac.var_names = [str(v) + '_' + adata_atac.var['gene_annotation'][i] for i,v in enumerate(adata_atac.var_names)]
adata_atac = scanpy_recipe(adata_atac,False,resolutions=[0.5,1,2,3,4,5,6],modality='atac',pca_n_comps=100,n_top_genes=60000)

adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_adt = sc.read('adata_after_scanpy_recipe_adt_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_atac = sc.read('adata_after_scanpy_recipe_atac_0.5_1_2_3_4_5_6_umap_True.h5ad')

# combine
adata_cite = concat_rna_and_other(adata_rna,adata_adt,'other','adt','AB_')
adata_tea = concat_rna_and_other(adata_cite,adata_atac,'rna','atac','')
adata_tea = make_sure_adata_writable(adata_tea,True)
adata_tea.write('input_adt_umap.h5ad')

'''
the input_adt_umap.h5ad combined AnnData as the input for scTriangulate, this file is on https://www.synapse.org/#!Synapse:syn26320350
'''

# plot
all_rs = ['sctri_' + m + '_leiden_' + str(i) for m in ['rna','atac','adt'] for i in [0.5,1,2,3,4,5,6]]
all_rs.extend(['azimuth','doublet_scores'])
umap_dual_view_save(adata_tea,all_rs)

# run
adata_tea = sc.read('input_adt_umap.h5ad')
all_rs = ['sctri_' + m + '_leiden_' + str(i) for m in ['rna','atac','adt'] for i in [1,2,3]]
sctri = ScTriangulate(dir='output_one',adata=adata_tea,query=all_rs,add_metrics={},predict_doublet='precomputed')
sctri.lazy_run(scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_adt_leiden_1','sctri_atac_leiden_1'])

'''
the result file after_pruned_assess.p will be on https://www.synapse.org/#!Synapse:syn26320350
'''

# insights
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
features = ['AB_CD56','CLEC10A','CLEC5A','CLEC4D','HLA-DRA','S100A9','FCGR3A']
peaks = [item for item in adata_atac.var_names if 'CASC15' in item]
sc.pl.umap(sctri.adata,color=features+peaks,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('output_one/previous_conclusions.pdf',bbox_inches='tight')
plt.close()

## modality contribution and confidence
sctri.modality_contributions()
for col in ['adt_contribution','atac_contribution','rna_contribution']:
    sctri.plot_umap(col,'continuous',umap_cmap='viridis')
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')


for cluster in ['sctri_adt_leiden_1@0','sctri_atac_leiden_1@2','sctri_rna_leiden_3@27','sctri_rna_leiden_1@2','sctri_rna_leiden_3@25']:
    sctri.plot_multi_modal_feature_rank(cluster)

for cluster in ['sctri_rna_leiden_2@12','sctri_rna_leiden_2@10','sctri_rna_leiden_1@1']:
    sctri.plot_multi_modal_feature_rank(cluster)

genes = ['CCR7','CCL5','TSHZ2','IFIT2','IL2RA']
sc.pl.umap(sctri.adata,color=genes,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('output_one/cd4_markers.pdf',bbox_inches='tight')
plt.close()

subset = ['sctri_adt_leiden_1@0','sctri_atac_leiden_1@2','sctri_rna_leiden_1@2']
sctri.plot_heterogeneity('azimuth','CD4 TCM','build',subset=subset)
for gene in ['ZBTB17','chr1_15948613_15951748_ZBTB17']:
    sctri.plot_heterogeneity('azimuth','CD4 TCM','single_gene',subset=subset,single_gene=gene,cmap='viridis')
for gene in ['STAM','ICOS','RUNX3','ROR2','AB_CD95']:
    sctri.plot_heterogeneity('azimuth','CD4 TCM','single_gene',subset=subset,single_gene=gene,cmap='viridis')
sys.exit('stop')

adts = []
for feature in sctri.adata.var_names:
    if feature.startswith('AB_'):
        adts.append(feature)
sc.pl.umap(sctri.adata,color=adts,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('output_one/all_adts.pdf')
plt.close()

subset = ['sctri_adt_leiden_2@20','sctri_rna_leiden_2@17']
for gene in ['HLA-DRA','S100A9','CLEC4D','CLEC5A','CLEC10A','FCGR3A','CD14','AB_CD16','AB_CD14']:
    sctri.plot_heterogeneity('azimuth','CD16 Mono','single_gene',subset=subset,single_gene=gene,cmap='viridis')




