#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os,sys
from sklearn.preprocessing import Binarizer
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import *


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# # inspect data
# adata_rna = sc.read_10x_mtx(path='./data/',prefix='GSM3309836_MF05.')

# # calr = pd.read_csv('GSM3309842_MF05.CALR.txt',sep='\t',index_col=0)['num.MUT.call']
# # nfe2 = pd.read_csv('GSM3309843_MF05.NFE2.txt',sep='\t',index_col=0)['num.MUT.call']
# # sf3b1 = pd.read_csv('GSM3309844_MF05.SF3B1.txt',sep='\t',index_col=0)['num.MUT.call']
# # amplicons = pd.concat((calr,nfe2,sf3b1),axis=1,join='outer',keys=('calr','nfe2','sf3b1'))
# # amplicon_annotation = []
# # for row in amplicons.itertuples():
# #     if row.nfe2 > 0:
# #         amplicon_annotation.append('NFE2_CALR_SF3B1_MUT')
# #     elif row.calr > 0:
# #         amplicon_annotation.append('CALR_SF3B1_MUT')
# #     elif row.sf3b1 > 0:
# #         amplicon_annotation.append('SF3B1_MUT')
# #     else:
# #         if row.sf3b1 == 0:
# #             amplicon_annotation.append('WT')
# #         elif np.isnan(row.sf3b1):   # row.sf3b1 is a float na
# #             amplicon_annotation.append('SF3B1_MUT')  # presumbly
# # amplicons['amplicon_annotation'] = amplicon_annotation
# # amplicons.to_csv('amplicon.txt',sep='\t')

# amplicon = pd.read_csv('amplicon_nathan.txt',sep='\t',index_col=0)
# umaps = pd.read_csv('exp.Centroids-HCA-CD34+__GSM3309836_MF05.matrix_matrix_CPTT-filtered-ReOrdered-UMAP_scores.txt',sep='\t',index_col=0,header=None)
# umaps.columns = ['umap_x','umap_y']
# common = list(set(adata_rna.obs_names).intersection(set(amplicon.index)).intersection(set(umaps.index)))
# adata_rna = adata_rna[common,:]
# amplicon = amplicon.loc[common,:]
# umaps = umaps.loc[common,:]
# add_umap(adata_rna,umaps,mode='pandas_memory',cols=['umap_x','umap_y'])

# # expand the amplicon df and construct another adata
# tmp = []
# for item in amplicon['Assigned Genotype']:
#     if item == 'NFE2-CALR-SF3B1-MUT':
#         tmp.append((1,1,1))
#     elif item == 'SF3B1-CALR-MUT':
#         tmp.append((0,1,1))
#     elif item == 'SF3B1-MUT':
#         tmp.append((0,0,1))
#     elif item == 'SF3B1-WT':
#         tmp.append((0,0,0))
# amplicon['anno_NFE2'] = list(zip(*tmp))[0]
# amplicon['anno_CALR'] = list(zip(*tmp))[1]
# amplicon['anno_SF3B1'] = list(zip(*tmp))[2]
# adata_amplicon = ad.AnnData(X=amplicon.loc[:,['anno_NFE2','anno_CALR','anno_SF3B1']].values,var=pd.DataFrame(index=['anno_NFE2','anno_CALR','anno_SF3B1']),obs=amplicon)

# # clustering on rna
# # adata_rna = scanpy_recipe(adata_rna,False,[0.5,1,2,3,4,5],'rna',False,pca_n_comps=50,n_top_genes=3000)
# adata_rna = sc.read('adata_after_scanpy_recipe_rna_0.5_1_2_3_4_5_umap_False.h5ad')
# # cols = ['sctri_rna_leiden_{}'.format(i) for i in [0.5,1,2,3,4,5]]
# # umap_dual_view_save(adata_rna,cols)

# run scTriangulate
# adata_combined = concat_rna_and_other(adata_rna,adata_amplicon,'rna','amplicon','amplicon_')
# sctri = ScTriangulate(dir='output_mut_two',adata=adata_combined,query=['Assigned Cell-Type','Genotype-Cell-Type'])
# sctri.lazy_run(viewer_heterogeneity_keys=['Assigned_Cell-Type'])
sctri = ScTriangulate.deserialize('output_mut_two/after_pruned_assess.p')
# sctri.adata.obs.drop(columns=['anno_NFE2','anno_CALR','anno_SF3B1'])
# sc.pp.highly_variable_genes(sctri.adata,flavor='seurat',n_top_genes=3000)
# hv_genes = sctri.adata.var['highly_variable'].loc[sctri.adata.var['highly_variable']].index.tolist()
# hv_features = list(set(hv_genes + ['amplicon_{}'.format(item) for item in ['anno_NFE2','anno_CALR','anno_SF3B1']]))
# adata_nca = nca_embedding(sctri.adata,10,'pruned','umap',hv_features=hv_features)
# adata_nca.write('output_mut_two/adata_nca.h5ad')
adata_nca = sc.read('output_mut_two/adata_nca.h5ad')
# sctri.adata.obs['pruned'].to_csv('output_mut_two/barcode2pruned.txt',sep='\t')
# pd.DataFrame(data=adata_nca.obsm['X_umap'],index=adata_nca.obs_names,columns=['umap_x','umap_y']).to_csv('output_mut_two/nca_umap_coords.txt',sep='\t')
# sc.pl.umap(adata_nca,color=['amplicon_{}'.format(item) for item in ['anno_NFE2','anno_CALR','anno_SF3B1']],color_map=bg_greyed_cmap('viridis'));plt.savefig('output_mut_two/nca_plot_mutation.pdf',bbox_inches='tight');plt.close()
# sctri.adata.obsm['X_umap'] = adata_nca.obsm['X_umap']
# sctri.modality_contributions(regex_dict={'mutation':r'amplicon_'})
# for col in ['confidence','rna_contribution','mutation_contribution']:
#     sctri.plot_umap(col,'continuous',umap_cmap=bg_greyed_cmap('viridis'))
# for col in ['pruned','Assigned_Cell-Type','final_annotation']:
#     sctri.plot_umap(col,'category')
# sctri.gene_to_df(mode='marker_genes',key='pruned')
# sctri.plot_multi_modal_feature_rank(cluster='Genotype-Cell-Type@ERP-Early__NFE2-CALR-SF3B1-MUT',regex_dict={'mutation':r'amplicon_'})
# sctri.plot_heterogeneity('Assigned_Cell-Type','ERP-Early','sankey')
# sctri.plot_two_column_sankey('Assigned_Cell-Type','pruned',0.5)
# sctri.plot_two_column_sankey('Assigned_Cell-Type','pruned',opacity=0.5,text=True,margin=5)
ScTriangulate.salvage_run('build_all_viewers','output_mut_two/after_pruned_assess.p',viewer_cluster=False,
                           viewer_heterogeneity_keys=['Assigned_Cell-Type'],other_umap=adata_nca.obsm['X_umap'],
                           heatmap_scale='mean',heatmap_cmap=retrieve_pretty_cmap('altanalyze'),heatmap_cbar_scale=0.15)
# plot_coexpression(sctri.adata,gene2='amplicon_anno_NFE2',gene1='NFE2',kind='contour',contour_train_kde='semi_valid',outdir='output_mut_two')













