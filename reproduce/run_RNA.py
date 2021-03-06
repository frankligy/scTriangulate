#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os
import sys
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap,generate_gradient
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


'''
All the inputs and intermediate files are deposited on https://www.synapse.org/#!Synapse:syn26320659
'''

# this dataset is downloaded from
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3?

# load the data
adata = sc.read_10x_h5('./pbmc_10k_v3_filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata.obs.to_csv('check_obs.txt',sep='\t')


# qc
for key in ['n_genes_by_counts','total_counts','pct_counts_mt']:
    sc.pl.violin(adata,key,jitter=0.4)
    plt.savefig('qc_violin_{}.pdf'.format(key),bbox_inches='tight')
    plt.close()

sc.pl.scatter(adata,x='n_genes_by_counts',y='total_counts',color='pct_counts_mt')
plt.savefig('qc_scatter.pdf',bbox_inches='tight')
plt.close()


# let's tentatively suggests 20% mitochondria, 500 min counts, 300 min genes
print(adata)   # 11769 × 33538
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_cells(adata, min_counts=500)
adata = adata[adata.obs.pct_counts_mt < 20, :]
print(adata)   # 11022 × 33538
adata.write('post_qc.h5ad')  # run Azimuth with this object

# get the clustering and umap
adata = sc.read('post_qc.h5ad')
doublet_map = doublet_predict(adata)
adding_azimuth(adata,'azimuth_pred.tsv')
adata = scanpy_recipe(adata,is_log=False,resolutions=[0.5,1,2,3,4,5,6],pca_n_comps=50,n_top_genes=3000)
adata.obs.to_csv('check_obs.txt',sep='\t')
cols = adata.obs.columns.tolist()
umap_dual_view_save(adata,cols)

# run scTriangulate
adata = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_6_umap_True.h5ad')
sctri = ScTriangulate(dir='./output_two',adata=adata,query=['sctri_rna_leiden_1','sctri_rna_leiden_2','sctri_rna_leiden_3'],predict_doublet='precomputed')
sctri.lazy_run(viewer_heterogeneity_keys=[sctri.reference,'azimuth'])

# insights
sctri = ScTriangulate.deserialize('output_two/after_pruned_assess.p')

'''overall insights'''
sctri.obs_to_df()

'''looking from leiden1'''
# cluster0, CD4 TCM and CD4 TEM, sctriangulate success
subset = ['sctri_rna_leiden_1@0','sctri_rna_leiden_3@2']
sctri.plot_heterogeneity('sctri_rna_leiden_1','0','violin',subset=subset,genes=['CCR7','CCL5'])
sctri.plot_heterogeneity('sctri_rna_leiden_1','0','dual_gene',subset=subset,dual_gene=['CCR7','CCL5'])
sctri.plot_heterogeneity('sctri_rna_leiden_1','0','umap',subset=subset)

# cluster1,4,5,12,17, CD14 Monocyte, azimuth undercluster, sctriangulate identify new populations
subset = ['sctri_rna_leiden_2@0','sctri_rna_leiden_3@10','sctri_rna_leiden_3@7','sctri_rna_leiden_1@4','sctri_rna_leiden_1@5','sctri_rna_leiden_3@13','sctri_rna_leiden_1@17']
merge = [('sctri_rna_leiden_1@5','sctri_rna_leiden_3@13'),('sctri_rna_leiden_2@0','sctri_rna_leiden_3@10','sctri_rna_leiden_3@7','sctri_rna_leiden_1@4')]
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CD14',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='FCGR3A',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC10A',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC5A',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='CLEC4D',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='MX1',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='MX2',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='IFI44',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='S100A9',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','single_gene',subset=subset,single_gene='HLA-DRA',cmap='viridis')
sctri.plot_heterogeneity('azimuth','CD14 Mono','umap',subset=subset)   
sctri.plot_heterogeneity('azimuth','CD14 Mono','build',subset=subset)
sys.exit('stop')
sctri.plot_heterogeneity('azimuth','CD14 Mono','umap',subset=subset,merge=merge)

# cluster2, CD4 naive, same

# cluster3, should be B memory + B intermediate, but none of the leiden recover it, I think it is due to the marker is not stable
sc.pl.umap(sctri.adata,color=['MS4A1','TNFRSF13B','IGHM','IGHD','AIM2','CD79A','LINC01857','RALGPS2','BANK1','CD79B'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('check_b.pdf',bbox_inches='tight')
plt.close()

sc.pl.umap(sctri.adata,color=['MS4A1','COCH','AIM2','BANK1','SSPN','CD79A','TEX9','RALGPS2','TNFRSF13C','LINC01781'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('check_b1.pdf',bbox_inches='tight')
plt.close()

# cluster6, MAIT + gammaT, success
sctri.plot_heterogeneity('sctri_rna_leiden_1','6','umap')
sctri.plot_heterogeneity('sctri_rna_leiden_1','6','dual_gene',dual_gene=['SLC4A10','TRDC'])

# cluster7, NK, cd56 fail
sctri.plot_heterogeneity('sctri_rna_leiden_1','7','umap')
sctri.plot_heterogeneity('sctri_rna_leiden_1','7','single_gene',single_gene='NCAM1',cmap='viridis')

# cluster8, CD8 TCM and TEM
subset = ['sctri_rna_leiden_2@10','sctri_rna_leiden_3@18']
sctri.plot_heterogeneity('sctri_rna_leiden_1','8','umap',subset=subset)
sctri.plot_heterogeneity('sctri_rna_leiden_1','8','dual_gene',subset=subset,dual_gene=['IL7R','GZMH'])

# cluster9, B naive, same

# cluster10, CD16 Mono, same

# cluster11, CD8 naive

# cluster12, potentially doublet, predicted by scrublet with high confidence, also no obvious marker genes

# cluster13 has three parts, Treg, proliferating, plasmablast
sctri.plot_heterogeneity('sctri_rna_leiden_1','13','umap')
sctri.plot_heterogeneity('sctri_rna_leiden_1','13','multi_gene',multi_gene=['MZB1','IL2RA','MKI67'])

# cluster 14, potentially doublet, also no obvious marker genes

# cluster 15, platelet, same

# cluster 16, since we will remove ASDC, so that's fine, same, just cDC2

# cluster 17, potentially doublet, no obvious marker genes

# cluster 19, pDC, same

# cluster 20, ILC + HSPC
sc.pl.umap(sctri.adata,color=['KIT','TRDC','TTLL10','LINC01229','SOX4','KLRB1','TNFRSF18','TNFRSF4', 'IL1R1', 'HPGDS'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
plt.savefig('check_ilc.pdf',bbox_inches='tight')
plt.close()

# cluster 21, cDC1, same


'''addtional plots for figure'''

# score_justify plot, using CD8 memory cell as an example
subset = ['sctri_rna_leiden_2@10','sctri_rna_leiden_3@18']
sctri.plot_heterogeneity('sctri_rna_leiden_1','8','umap',subset=subset)
subset = ['sctri_rna_leiden_3@18']
sctri.plot_heterogeneity('sctri_rna_leiden_3','18','umap',subset=subset)


# retrieve the scores from scTriangulate viewer, basically the metrics and shapley values for each cluster, and then plot
y1 = [0.69,0.41,0.86,0.49,0.02]
y2 = [0.90,0.58,0.93,0.66,6.67]
x1 = [(lambda x:3*x+1)(i) for i in range(5)]
x2 = [(lambda x:3*x+2)(i) for i in range(5)]
x = x1 + x2
y = y1 + y2
fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':[0.3,0.7],'hspace':0.1})
for ax in axes:
    ax.bar(x=x,height=y,color=np.repeat(['#1f77b4','#ff7f0e'],5),edgecolor='k')
    for i in range(len(x)):
        ax.text(x[i],y[i]+0.1,str(y[i]),va='center',ha='center')

axes[0].set_ylim([6,7])
axes[1].set_ylim([0,1])

axes[0].spines['bottom'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[0].tick_params(bottom=False)

d = 0.015
axes[0].plot((-d,d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[0].plot((1-d,1+d),(-d,d),transform=axes[0].transAxes,clip_on=False,color='k')
axes[1].plot((-d,d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')
axes[1].plot((1-d,1+d),(1-d,1+d),transform=axes[1].transAxes,clip_on=False,color='k')

t = [(lambda x:3*x+1.5)(i) for i in range(5)]
axes[1].set_xticks(t)
axes[1].set_xticklabels(['Reassign Score','SCCAF Score','TF-IDF10 Score','TF-IDF 5 Score','Shapley Value'],rotation=60)

axes[1].legend(handles=[mpatch.Patch(color=i) for i in ['#1f77b4','#ff7f0e']],labels=['leiden1@9','leiden3@7'],
          loc='upper left',bbox_to_anchor=(1,1),frameon=False)
axes[1].set_xlabel('Metrics and Shapley')
axes[1].set_ylabel('Value')

plt.savefig('score_justify.pdf',bbox_inches='tight')
plt.close()

# confidence plot
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')

# pruned_final_annotation
sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
sctri.plot_umap('pruned_final_annotation','category')

# plot concordance
confusion_df = sctri.plot_concordance(key1='azimuth',key2='pruned',style='3dbar')
confusion_df.to_csv('output_two/azimuth_pruned_confusion.txt',sep='\t')
generate_gradient(cmap='jet',name='jet')


# v-measure plot
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure
result = sctri.adata.obs
azimuth = LabelEncoder().fit_transform(result['azimuth'].values)
leiden1 = LabelEncoder().fit_transform(result['leiden1'].values)
leiden2 = LabelEncoder().fit_transform(result['leiden2'].values)
leiden3 = LabelEncoder().fit_transform(result['leiden3'].values)
pruned = LabelEncoder().fit_transform(result['pruned'].values)

for test in [leiden1,leiden2,leiden3,pruned]:
    print('ARI:{}'.format(adjusted_rand_score(azimuth,test)))
    print('V:{}'.format(homogeneity_completeness_v_measure(azimuth, test)))
    print('NMI:{}'.format(normalized_mutual_info_score(azimuth,test)))

'''
ARI:0.5011579040293368
V:(0.8117030890841928, 0.7045263914413321, 0.7543267765458033)
NMI:0.7543267765458033
ARI:0.4328692303679871
V:(0.8248995749634627, 0.6655871491412404, 0.7367292140861539)
NMI:0.736729214086154
ARI:0.30173849414833653
V:(0.8381897816150324, 0.5934047555759651, 0.694869656965142)
NMI:0.6948696569651421
ARI:0.35295807833354215
V:(0.838036299746353, 0.6312244083998317, 0.7200750208482997)
NMI:0.7200750208482996
'''

fig,ax = plt.subplots()
ax.plot([1,2,3,4],[0.8117,0.8249,0.8382,0.8380],label='Homogeneity',marker='o',linestyle='--')
ax.plot([1,2,3,4],[0.7045,0.6655,0.5934,0.6312],label='Completeness', marker='o', linestyle='--')
ax.plot([1,2,3,4], [0.7543,0.7367,0.6948,0.7201], label='NMI', marker='o', linestyle='--')
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1,1))
ax.set_xticks([1,2,3,4])
ax.set_xticklabels(['r=1','r=2','r=3','scTriangulate'])
plt.savefig('check.pdf',bbox_inches='tight')
plt.close()




'''for building the tutorial'''
# sctri = ScTriangulate.deserialize('output_one/after_pruned_assess.p')
# subset = ['sctri_rna_leiden_2@7','sctri_rna_leiden_2@18']
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='umap',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='build',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='violin',genes=['TRDC','SLC4A10'],subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='single_gene',single_gene='TRDC',subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='dual_gene',dual_gene=['TRDC','SLC4A10'],subset=subset)
# sctri.plot_heterogeneity(key='sctri_rna_leiden_1',cluster='6',style='sankey',subset=subset)

# sctri.plot_confusion(name='confusion_reassign',key='sctri_rna_leiden_1')
# sctri.plot_circular_barplot(key='sctri_rna_leiden_1',col='pruned')

# sctri.get_metrics_and_shapley(barcode='AAACCCACATCCAATG-1',save=True)






















