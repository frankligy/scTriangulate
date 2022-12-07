#! /data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,sys
import scanpy as sc
import matplotlib.pyplot as plt

from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap

'''
All the input, intermediate files are deposited on https://www.synapse.org/#!Synapse:syn26320567
'''

# load the data
adata = sc.read_10x_h5('28WM_ND19-341__TNC-RNA-ADT.h5',gex_only=False)
adata_rna = adata[:,adata.var['feature_types']=='Gene Expression']
adata_adt = adata[:,adata.var['feature_types']=='Antibody Capture']  # 8491

adata_rna.var_names_make_unique()
adata_adt.var_names_make_unique()

# qc, rna
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
    sc.pl.violin(adata_rna,key,jitter=0.4)
    plt.savefig('qc_rna_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_rna,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
plt.savefig('qc_rna_scatter.pdf',bbox_inches='tight')
plt.close()

# qc, adt
sc.pp.calculate_qc_metrics(adata_adt, var_type='adts',percent_top=None, log1p=False, inplace=True)

for key in ['n_adts_by_counts','total_counts']:
    sc.pl.violin(adata_adt,key,jitter=0.4)
    plt.savefig('qc_adt_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata_adt,x='n_adts_by_counts',y='total_counts')
plt.savefig('qc_adt_scatter.pdf',bbox_inches='tight')
plt.close()

# decision, min_genes = 300, min_counts= 500, mito < 20
sc.pp.filter_cells(adata_rna, min_genes=300)
sc.pp.filter_cells(adata_rna, min_counts=500)
adata_rna = adata_rna[adata_rna.obs.pct_counts_mt < 20, :]
adata_adt = adata_adt[adata_rna.obs_names,:]   # 6406

adata_rna.write('post_rna_qc.h5ad')

# get clusters and input
doublet_map = doublet_predict(adata_rna)
adding_azimuth(adata_rna,'azimuth_pred.tsv')
adata_rna = scanpy_recipe(adata_rna,False,resolutions=[0.5,1,2,3,4,5,6],modality='rna',pca_n_comps=50)
adata_adt = scanpy_recipe(adata_adt,False,resolutions=[0.5,1,2,3,4,5,6],modality='adt',pca_n_comps=15)

adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_adt = sc.read('adata_after_scanpy_recipe_adt_0.5_1_2_3_4_5_6_umap_True.h5ad')
adata_combine = concat_rna_and_other(adata_rna,adata_adt,umap='other',name='adt',prefix='AB_')

# plot
all_r = ['sctri_{}_leiden_{}'.format(m,i) for i in [0.5,1,2,3,4,5,6] for m in ['rna','adt']]
cols = all_r + ['azimuth','doublet_scores']
umap_dual_view_save(adata_combine,cols)

adata_combine = make_sure_adata_writable(adata_combine,delete=True)
adata_combine.write('combined_rna_adt.h5ad')



# run scTriangulate
adata_combine = sc.read('combined_rna_adt.h5ad')
sctri = ScTriangulate(dir='output_one',adata=adata_combine,add_metrics={},predict_doublet='precomputed',
                      query=['sctri_adt_leiden_1','sctri_adt_leiden_2','sctri_adt_leiden_3','sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'])
sctri.lazy_run(viewer_heterogeneity_keys=['azimuth','sctri_adt_leiden_1','sctri_rna_leiden_1'])



# insights on TNC1
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
sc.pl.umap(sctri.adata,color=['AB_CD2','AB_CD7','AB_CD14','AB_CD49f','AB_CD71'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('correspond_adts.pdf',bbox_inches='tight')
plt.close()
sctri.plot_umap('AB_CD56','continuous',umap_cmap='viridis')

## WNN
wnn = pd.read_csv('wnn_metadata.txt',sep='\t',index_col=0)
for col in ['wsnn_res.0.25','wsnn_res.1','wsnn_res.2','wsnn_res.3','wsnn_res.3.75']:
    mapping = wnn[col].to_dict()
    sctri.adata.obs[col] = [str(item) for item in sctri.adata.obs_names.map(mapping).values]
    sctri.plot_umap(col,'category')

## mix-and-match schema
col = []
for item in sctri.adata.obs['pruned']:
    if 'leiden_1@' in item:
        col.append('resolution1')
    elif 'leiden_2@' in item:
        col.append('resolution2')
    elif 'leiden_3@' in item:
        col.append('resolution3')
sctri.adata.obs['resolution_distribution'] = col
sctri.plot_umap('resolution_distribution','category')


## plot all adts
all_adts = ['AB_' + item for item in adata_adt.var_names.tolist()]
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
sc.pl.umap(sctri.adata,color=all_adts,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('all_adt.pdf',bbox_inches='tight')
plt.close()

## output the obs
sctri.obs_to_df()

## adt and rna contributions
sctri.modality_contributions()
for col in ['adt_contribution','rna_contribution']:
    sctri.plot_umap(col,'continuous',umap_cmap='viridis')


## CD8, plot rank
for cluster in ['sctri_adt_leiden_3@32','sctri_rna_leiden_1@13','sctri_rna_leiden_3@6','sctri_adt_leiden_3@37','sctri_rna_leiden_1@9']:
    sctri.plot_multi_modal_feature_rank(cluster=cluster)

# Mono population
genes = ['HLA-DRA','S100A9','CLEC10A','CLEC5A','CLEC4D','FCGR3A','MX1','MX2','IFI44']
sc.pl.umap(sctri.adata,color=genes,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('mono_markers.pdf',bbox_inches='tight')
plt.close()

subset = ['sctri_rna_leiden_1@4','sctri_rna_leiden_3@15','sctri_rna_leiden_1@0','sctri_rna_leiden_1@11','sctri_rna_leiden_3@11']
sctri.plot_heterogeneity('azimuth','CD14 Mono','build',subset=subset)

# confidence
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')


# metrics

from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure
result = sctri.adata.obs
azimuth = LabelEncoder().fit_transform(result['azimuth'].values)
rna_leiden1 = LabelEncoder().fit_transform(result['sctri_rna_leiden_1'].values)
rna_leiden2 = LabelEncoder().fit_transform(result['sctri_rna_leiden_2'].values)
rna_leiden3 = LabelEncoder().fit_transform(result['sctri_rna_leiden_3'].values)
adt_leiden1 = LabelEncoder().fit_transform(result['sctri_adt_leiden_1'].values)
adt_leiden2 = LabelEncoder().fit_transform(result['sctri_adt_leiden_2'].values)
adt_leiden3 = LabelEncoder().fit_transform(result['sctri_adt_leiden_3'].values)
pruned = LabelEncoder().fit_transform(result['pruned'].values)

for test in [rna_leiden1,rna_leiden2,rna_leiden3,adt_leiden1,adt_leiden2,adt_leiden3,pruned]:
    print('ARI:{}'.format(adjusted_rand_score(azimuth,test)))
    print('V:{}'.format(homogeneity_completeness_v_measure(azimuth, test)))
    print('NMI:{}'.format(normalized_mutual_info_score(azimuth,test)))

'''
ARI:0.5596547814072221
V:(0.8015564188680933, 0.7127135838054284, 0.7545287787582127)
NMI:0.7545287787582127

ARI:0.49240529240720754
V:(0.8222715464906468, 0.6726673752279325, 0.7399837342740611)
NMI:0.739983734274061

ARI:0.33244409885587145
V:(0.8362610940736923, 0.5802281860552212, 0.6851054427164793)
NMI:0.6851054427164793

ARI:0.6333414470598028
V:(0.7248657478268431, 0.6718498810771979, 0.6973516389396073)
NMI:0.6973516389396073

ARI:0.39297006413473284
V:(0.7600388434044909, 0.5827727214256064, 0.6597052285829074)
NMI:0.6597052285829075

ARI:0.29778256926907254
V:(0.780396540691627, 0.5419996009053327, 0.6397093885830141)
NMI:0.6397093885830142

ARI:0.5866793579034331
V:(0.8211540548316435, 0.7201708388194243, 0.7673543805125276)
NMI:0.7673543805125275
'''

fig,ax = plt.subplots()
ax.plot([1,2,3,4,5,6,7],[0.8015,0.8222,0.8362,0.7248,0.7600,0.7803,0.8215],label='Homogeneity',marker='o',linestyle='--')
ax.plot([1,2,3,4,5,6,7],[0.7127,0.6726,0.5802,0.6718,0.5827,0.5419,0.7201],label='Completeness', marker='o', linestyle='--')
ax.plot([1,2,3,4,5,6,7], [0.7545,0.7399,0.6851,0.6973,0.6597,0.6397,0.7673], label='NMI', marker='o', linestyle='--')
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
ax.set_xticks([1,2,3,4,5,6,7])
ax.set_xticklabels(['rna1','rna2','rna3','adt1','adt2','adt3','scTriangulate'])
plt.savefig('check.pdf',bbox_inches='tight')
plt.close()








