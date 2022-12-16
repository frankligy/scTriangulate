#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
import os,sys
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


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
adata_rna.write('post_rna_qc.h5ad')
adata_adt.write('post_adt_qc.h5ad')
adata_atac.write('post_atac_qc.h5ad')



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
adata_tea = concat_rna_and_other(adata_cite,adata_atac,'other','atac','')
adata_tea = make_sure_adata_writable(adata_tea,True)
adata_tea.write('input_atac_umap.h5ad')

# plot
all_rs = ['sctri_' + m + '_leiden_' + str(i) for m in ['rna','atac','adt'] for i in [0.5,1,2,3,4,5,6]]
all_rs.extend(['azimuth','doublet_scores'])
umap_dual_view_save(adata_tea,all_rs)

# run
adata_tea = sc.read('input_adt_umap.h5ad')
all_rs = ['sctri_' + m + '_leiden_' + str(i) for m in ['rna','atac','adt'] for i in [1,2,3]]
sctri = ScTriangulate(dir='output_one',adata=adata_tea,query=all_rs,add_metrics={},predict_doublet='precomputed')
sctri.lazy_run(scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_adt_leiden_1','sctri_atac_leiden_1'])

sctri = ScTriangulate(dir='output_two',adata=adata_tea,query=all_rs,predict_doublet='precomputed')
sctri.lazy_run(scale_sccaf=False,viewer_heterogeneity_keys=['azimuth','sctri_rna_leiden_1','sctri_adt_leiden_1','sctri_atac_leiden_1'])

# salvage run
ScTriangulate.salvage_run(step_to_start='assess_pruned',last_step_file='output_one/after_rank_pruning.p',win_fraction_cutoff=0.4,viewer_cluster=False,
                          viewer_heterogeneity_keys=['azimuth'])



# insights
sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
sctri.adata.obs.to_csv('to_bm_ensemble.txt',sep='\t')

# for ensemble
cols = ['hgpa','mcla','hbgf']
add_annotations(sctri.adata,'./from_bm_ensemble.txt',cols,0,cols)
sctri.cluster_performance(cluster='pruned',competitors=cols,reference='azimuth',show_cluster_number=True,metrics=True)
for col in cols:
    sctri.adata.obs[col] = sctri.adata.obs[col].astype('str')
    sctri.adata.obs[col] = sctri.adata.obs[col].astype('category')
    sctri.plot_umap(col=col,kind='category')




# for wnn or individual
cols = ['wsnn_res.0.25','wsnn_res.0.5','wsnn_res.0.75','wsnn_res.1','wsnn_res.1.25','wsnn_res.1.5','wsnn_res.1.75',
        'wsnn_res.2','wsnn_res.2.25','wsnn_res.2.5','wsnn_res.2.75','wsnn_res.3','wsnn_res.3.25','wsnn_res.3.5','wsnn_res.3.75']
add_annotations(sctri.adata,'./wnn_metadata.txt',cols,0,cols)
sctri.cluster_performance(cluster='pruned',competitors=cols,reference='azimuth',show_cluster_number=True,metrics=True)
sctri.cluster_performance(cluster='pruned',competitors=sctri.query,reference='azimuth',show_cluster_number=True,metrics=True)


# for totalVI
cols = ['leiden_totalVI_{}'.format(r) for r in [0.5,1,2,3,4,5,6]]
add_annotations(sctri.adata,'./totalVI_obs.txt',cols,0,cols)
sctri.cluster_performance(cluster='pruned',competitors=cols,reference='azimuth',show_cluster_number=True,metrics=True)



features = ['AB_CD56','CLEC10A','CLEC5A','CLEC4D','HLA-DRA','S100A9','FCGR3A']
peaks = [item for item in adata_atac.var_names if 'CASC15' in item]
sc.pl.umap(sctri.adata,color=features+peaks,cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('output_one/previous_conclusions.pdf',bbox_inches='tight')
plt.close()

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
# sctri.plot_heterogeneity('azimuth','CD4 TCM','build',subset=subset)
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

# MAIT, NKT
sctri.plot_cluster_feature('azimuth','MAIT','location')



'''
answer reviewer questions (why shapley versus rank, why 0.01, why pruning and reassigning)
original, shapley_all_or_none
rank
rank_all_or_none
shapley
'''
for opt in ['rank_all_or_none','rank','shapley']:
    os.mkdir('output_one_{}'.format(opt))
    subprocess.run(['cp','output_one/after_metrics.p','output_one_{}'.format(opt)])
    ScTriangulate.salvage_run(step_to_start='run_shapley',last_step_file='output_one/after_metrics.p',outdir='output_one_{}'.format(opt),
                              shapley_mode=opt,shapley_bonus=0.01,assess_pruned=False,viewer_cluster=False,viewer_heterogeneity=False)

sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
for opt in ['rank_all_or_none','rank','shapley']:
    new_sctri = ScTriangulate.deserialize('output_one_{}/after_rank_pruning.p'.format(opt))
    new_sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(new_sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
    new_sctri.add_to_invalid_by_win_fraction(percent=0.25)
    new_sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=new_sctri.reference)
    new_obs = new_sctri.adata.obs
    add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_{}'.format(opt)],'\t','memory')

sctri.cluster_performance(cluster='pruned',competitors=['pruned_{}'.format(opt) for opt in ['rank_all_or_none','rank','shapley']],
                          reference='azimuth',show_cluster_number=True,metrics=True)

for bonus in [0,0.005,0.01,0.05,0.1]:
    os.mkdir('output_one_bonus_{}'.format(bonus))
    subprocess.run(['cp','output_one/after_metrics.p','output_one_bonus_{}'.format(bonus)])
    ScTriangulate.salvage_run(step_to_start='run_shapley',last_step_file='output_one/after_metrics.p',outdir='output_one_bonus_{}'.format(bonus),
                              shapley_mode='shapley_all_or_none',shapley_bonus=bonus,assess_pruned=False,viewer_cluster=False,viewer_heterogeneity=False)

sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
for bonus in [0,0.005,0.01,0.05,0.1]:
    new_sctri = ScTriangulate.deserialize('output_one_bonus_{}/after_rank_pruning.p'.format(bonus))
    new_sctri.uns['raw_cluster_goodness'].to_csv(os.path.join(new_sctri.dir,'raw_cluster_goodness.txt'),sep='\t')
    new_sctri.add_to_invalid_by_win_fraction(percent=0.25)
    new_sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=new_sctri.reference)
    new_obs = new_sctri.adata.obs
    add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_bonus_{}'.format(bonus)],'\t','memory')

sctri.cluster_performance(cluster='pruned',competitors=['pruned_bonus_{}'.format(bonus) for bonus in [0,0.005,0.01,0.05,0.1]],
                          reference='azimuth',show_cluster_number=True,metrics=True)



sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
for wfc in [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]:
    sctri.clear_invalid()
    sctri.add_to_invalid_by_win_fraction(percent=wfc)
    sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=sctri.reference)
    new_obs = sctri.adata.obs
    add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_wfc_{}'.format(wfc)],'\t','memory')

sctri.adata.obs['pruned'] = sctri.adata.obs['pruned_wfc_0.25']
sctri.cluster_performance(cluster='pruned',competitors=['pruned_wfc_{}'.format(wfc) for wfc in [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6]],
                          reference='azimuth',show_cluster_number=True,metrics=True)

sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
for nabs in [1,5,10,15,20,25,30,35,40,45,50]:
    sctri.clear_invalid()
    sctri.add_to_invalid_by_win_fraction(percent=0.25)
    sctri.pruning(method='reassign',abs_thresh=nabs,remove1=True,reference=sctri.reference)
    new_obs = sctri.adata.obs
    add_annotations(sctri.adata,new_obs,['pruned'],0,['pruned_nabs_{}'.format(nabs)],'\t','memory')

sctri.adata.obs['pruned'] = sctri.adata.obs['pruned_nabs_10']
sctri.cluster_performance(cluster='pruned',competitors=['pruned_nabs_{}'.format(nabs) for nabs in [1,5,10,15,20,25,30,35,40,45,50]],
                          reference='azimuth',show_cluster_number=True,metrics=True)


