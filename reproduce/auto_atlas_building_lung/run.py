#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

'''
old env: #!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6
new env: #!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

old codebase: /data/salomonis2/LabFiles/Frank-Li/scTriangulate/src_backup
new codebase: /data/salomonis2/software
'''

import os
import sys
import scanpy as sc
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


'''
the mtx folder /data is on synapse
all the annotation txt file /group are on synapse
umap coordinate txt file is on synapse
assembLed input.h5ad (count matrix with umap coordinates and cluster annotations saved in obs) file is on synapse

the codes that are commented out have been used over the course of my development, I put them here just for transparency, 
you can ignore them to better follow the logics presented in the paper.
'''


# large_txt_to_mtx(int_file='GSE161382_counts_matrix.txt',out_folder='./data',gene_is_index=True,type_convert_to='int16')
# adata = mtx_to_adata(int_folder='./data',gene_is_index=True,feature='genes')
# adata = just_log_norm(adata)
# add_annotations(adata,'./Groups/groups.Columbia-round-2.txt',['label'],0,['Columbia'])
# add_annotations(adata,'./Groups/groups.Kaminsky-round-2.txt',['label'],0,['Kaminsky'])
# add_annotations(adata,'./Groups/groups.Krasnow-round-2.txt',['label'],0,['Krasnow'])
# add_annotations(adata,'./Groups/groups.Rhesus-round-2.txt',['label'],0,['Rhesus'])
# add_annotations(adata,'./Groups/groups.Sun.txt',['label'],0,['Sun'])
# add_umap(adata,'GSE161382_UMAP_coord.txt','pandas',['umap_x','umap_y'],0)
# adata.write('input.h5ad')
# umap_dual_view_save(adata,cols=['Columbia','Kaminsky','Krasnow','Rhesus','Sun'])

'''
Wang et al. (https://elifesciences.org/articles/62522) GSE161383 Sun dataset
Adam et al. (science.org/doi/10.1126/sciadv.aba1983) GSE136831 Kaminsky dataset

'''

adata = sc.read('input.h5ad')

'''opt1, with Rhesus, one tfidf'''
# sctri = ScTriangulate(dir='output_opt1_RY_one',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow','Rhesus'],add_metrics={})
# sctri.compute_metrics(parallel=True,scale_sccaf=True)
# sctri.serialize('after_metrics.p')
# sctri.compute_shapley(parallel=True)
# sctri.serialize('after_shapley.p')
# sctri.pruning(method='rank',scale_sccaf=True)
# sctri.serialize('after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.plot_umap('final_annotation','category')
# sctri.plot_umap('pruned','category')
# sctri.viewer_heterogeneity_figure(key='Sun')
# sctri.viewer_heterogeneity_html(key='Sun')
# sctri.viewer_cluster_feature_html()
# sctri.viewer_cluster_feature_figure()


# sctri = ScTriangulate.deserialize('output_opt1_RY_one/after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.plot_heterogeneity('pruned','Krasnow@Serous','umap')
# sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
# sctri.plot_umap('pruned_final_annotation','category')
# for col in ['Sun','Krasnow','Kaminsky','Columbia','Rhesus','final_annotation','pruned']:    
#     sctri.plot_umap(col,'category',format='png')
# hightlight1: AT1
# subset = ['Columbia@AT1','Kaminsky@ATI']
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','umap',subset=subset,format='png')
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','violin',subset=subset,format='pdf',genes=['IGF2BP2'])
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','single_gene',subset=subset,format='pdf',single_gene='IGF2BP2')
# hightlight2: Alveolar_macrophage
# subset = ['Kaminsky@Macrophage_Alveolar','Sun@Alveolar_macrophages']
# marker_dict = sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap',subset=subset,format='pdf')
# marker_dict['Kaminsky@Macrophage_Alveolar'].to_csv('Kaminsky.txt',sep='\t')
# marker_dict['Sun@Alveolar_macrophages'].to_csv('Sun.txt',sep='\t')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','umap',subset=subset,format='png')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','violin',subset=subset,format='pdf',genes=['INHBA','FABP4'])
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='INHBA')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='FABP4')
# marker_gene_dict = {'Kaminsky@Macrophage_Alveolar':['C1QB','C1QA','C1QC','INHBA','FABP4','HLA-DPA1','CD52','HLA-DPA1','HLA-DRB1'],
#                     'Sun@Alveolar_macrophages':['DOCK4','HIF1A','VCAN','ADAM9','STARD13','PLPP3','MERTK','GNA12']}
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)
# highlight3: dendritic cells
# subset = ['Kaminsky@cDC1','Krasnow@Plasmacytoid_Dendritic']
# sctri.plot_heterogeneity('Sun','dendritic_cells','umap',format='png',subset=subset)
# sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap',subset=subset)
# marker_gene_dict = {'Kaminsky@cDC1':['ASAP1','ATP2B1','CCSER1','SLC9A9','BCL6','CLNK','ID2','ZNF366','CFLAR','ENOX1','DAPP1','TOX','FMNL2',
#                     'SLCO3A1','RAB8B','CD83','FGD4','BASP1','ZBTB46','CPVL','SLC8A1','GPRIN3','SLCO5A1','CADM1','CCDC26'],
#                     'Krasnow@Plasmacytoid_Dendritic':['ZFAT','CCDC50','RUNX2','COBLL1','PTPRS','FAM129C','PLAC8',
#                     'LINC01478','PLXNA4','CUX2','ARHGAP24','NPC1','COL26A1','EPHB1','CLEC4C','COL24A1']}
# sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)
# sctri.plot_heterogeneity('Sun','dendritic_cells','violin',subset=subset,genes=['CADM1'])


'''opt1, with Rhesus, two tfidf'''
# sctri = ScTriangulate(dir='output_opt1_RY_two',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow','Rhesus'])
# sctri.compute_metrics(parallel=True,scale_sccaf=True)
# sctri.serialize('after_metrics.p')
# sctri.compute_shapley(parallel=True)
# sctri.serialize('after_shapley.p')
# sctri.pruning(method='rank',scale_sccaf=True)
# sctri.serialize('after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.plot_umap('final_annotation','category')
# sctri.plot_umap('pruned','category')
# sctri.viewer_heterogeneity_figure(key='Sun')
# sctri.viewer_heterogeneity_html(key='Sun')
# sctri.viewer_cluster_feature_html()
# sctri.viewer_cluster_feature_figure()

# sctri = ScTriangulate.deserialize('output_opt1_RY_two/after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
# sctri.plot_umap('pruned_final_annotation','category')

'''opt1, no Rhesus, one tfidf'''
# sctri = ScTriangulate(dir='output_opt1_RN_one',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow'],add_metrics={})
# sctri.compute_metrics(parallel=True,scale_sccaf=True)
# sctri.serialize('after_metrics.p')
# sctri.compute_shapley(parallel=True)
# sctri.serialize('after_shapley.p')
# sctri.pruning(method='rank',scale_sccaf=True)
# sctri.serialize('after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.plot_umap('final_annotation','category')
# sctri.plot_umap('pruned','category')
# sctri.viewer_heterogeneity_figure(key='Sun')
# sctri.viewer_heterogeneity_html(key='Sun')
# sctri.viewer_cluster_feature_html()
# sctri.viewer_cluster_feature_figure()

# sctri = ScTriangulate.deserialize('output_opt1_RN_one/after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
# sctri.plot_umap('pruned_final_annotation','category')

'''opt1, no Rhesus, two tfidf'''
### reproduce using new version
# sctri = ScTriangulate(dir='output_opt1_RN_two_new',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow'])
# sctri.lazy_run(win_fraction_cutoff=0.4)

### reproduce using latest version
# sctri = ScTriangulate(dir='output_opt1_RN_two_latest',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow'])
# sctri.lazy_run(win_fraction_cutoff=0.4,compute_metrics_parallel=False)  

sctri = ScTriangulate.deserialize('output_opt1_RN_two_latest/after_pruned_assess.p')
# sctri.extract_stability()
# # sctri.adata.obs.to_csv('output_opt1_RN_two_new/to_bm_ensemble.txt',sep='\t')
# cols = ['hgpa','mcla','hbgf']
# add_annotations(sctri.adata,'output_opt1_RN_two_new/from_bm_ensemble.txt',cols,0,cols)
# for col in cols:
#     sctri.adata.obs[col] = sctri.adata.obs[col].astype('str')
#     sctri.adata.obs[col] = sctri.adata.obs[col].astype('category')
#     sctri.plot_umap(col=col,kind='category')
# sys.exit('stop')

# add_annotations(sctri.adata,'azimuth_pred_lungmap.tsv',['predicted.celltype_level3'],0,['lungmap_azimuth_pred'],'\t','disk')
# sctri.adata.obs.loc[:,['Sun','Columbia','Kaminsky','Krasnow','pruned','lungmap_azimuth_pred']].to_csv('output_opt1_RN_two_latest/reviewer_compare_lungmap.txt',sep='\t')
# sctri.plot_umap(col='lungmap_azimuth_pred',kind='category')

# add_annotations(sctri.adata,'azimuth_pred_HCA.tsv',['predicted.ann_finest_level'],0,['HCA_azimuth_pred'],'\t','disk')
# sctri.adata.obs.loc[:,['Sun','Columbia','Kaminsky','Krasnow','pruned','HCA_azimuth_pred']].to_csv('output_opt1_RN_two_latest/reviewer_compare_HCA.txt',sep='\t')
# sctri.plot_umap(col='HCA_azimuth_pred',kind='category')


# sctri = ScTriangulate(dir='output_opt1_RN_two',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow'])
# sctri.compute_metrics(parallel=True,scale_sccaf=True)
# sctri.serialize('after_metrics.p')
# sctri.compute_shapley(parallel=True)
# sctri.serialize('after_shapley.p')
# sctri.pruning(method='rank',scale_sccaf=True)
# sctri.serialize('after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.plot_umap('final_annotation','category')
# sctri.plot_umap('pruned','category')
# sctri.viewer_heterogeneity_figure(key='Sun')
# sctri.viewer_heterogeneity_html(key='Sun')
# sctri.viewer_cluster_feature_html()
# sctri.viewer_cluster_feature_figure()

# sctri = ScTriangulate.deserialize('output_opt1_RN_two/after_prune_rank.p')
# sctri.add_to_invalid_by_win_fraction(percent=0.4)
# sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference='Sun')
# sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
# sctri.plot_umap('pruned_final_annotation','category')
# for col in ['Sun','Krasnow','Kaminsky','Columbia','Rhesus','final_annotation','pruned']:    
#     sctri.plot_umap(col,'category',format='png')

# sctri.plot_umap('confidence','continuous',umap_cmap='viridis')
# sys.exit('stop')

# hightlight1: AT1
# subset = ['Columbia@AT1','Kaminsky@ATI']
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','umap',subset=subset,format='png')
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','violin',subset=subset,format='pdf',genes=['IGFBP2'])
# marker_dict = sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','heatmap',subset=subset)
# marker_dict['Kaminsky@ATI'].to_csv('output_opt1_RN_two/Kaminsky@ATI.txt',sep='\t')
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','violin',subset=subset,format='pdf',genes=['IGF2BP2'])
# sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','single_gene',subset=subset,format='pdf',single_gene='IGF2BP2')

# # hightlight2: Alveolar_macrophage
# subset = ['Kaminsky@Macrophage_Alveolar','Sun@Alveolar_macrophages']
# marker_dict = sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap',subset=subset,format='pdf')
# marker_dict['Kaminsky@Macrophage_Alveolar'].to_csv('output_opt1_RN_two/Kaminsky.txt',sep='\t')
# marker_dict['Sun@Alveolar_macrophages'].to_csv('output_opt1_RN_two/Sun.txt',sep='\t')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','umap',subset=subset,format='png')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','violin',subset=subset,format='pdf',genes=['INHBA','FABP4'])
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='INHBA')
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='FABP4')
# marker_gene_dict = {'Kaminsky@Macrophage_Alveolar':['C1QB','C1QA','C1QC','INHBA','FABP4','HLA-DPA1','CD52','HLA-DPA1','HLA-DRB1'],
#                     'Sun@Alveolar_macrophages':['DOCK4','HIF1A','VCAN','ADAM9','STARD13','PLPP3','MERTK','GNA12']}
# sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)



# highlight3: dendritic cells
# sctri.adata.var_names.to_series().to_csv('output_opt1_RN_two_latest/all_genes.txt',sep='\t')
# subset = ['Kaminsky@cDC1','Krasnow@Plasmacytoid_Dendritic','Kaminsky@DC_Mature','Kaminsky@DC_Langerhans']
# marker_dict = sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap',subset=subset,format='pdf')
# marker_dict['Kaminsky@cDC1'].to_csv('output_opt1_RN_two_latest/Kaminsky@cDC1.txt',sep='\t')
# marker_dict['Kaminsky@DC_Mature'].to_csv('output_opt1_RN_two_latest/Kaminsky@DC_Mature.txt',sep='\t')
# marker_dict['Kaminsky@DC_Langerhans'].to_csv('output_opt1_RN_two_latest/Kaminsky@DC_Langerhans.txt',sep='\t')
# marker_dict['Krasnow@Plasmacytoid_Dendritic'].to_csv('output_opt1_RN_two_latest/Krasnow@Plasmacytoid_Dendritic.txt',sep='\t')
# sctri.plot_heterogeneity('Sun','dendritic_cells','umap',format='png',subset=subset)
# sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap',subset=subset)
# marker_gene_dict = {'Krasnow@Plasmacytoid_Dendritic':['TCF4','ZFAT','CCDC50','RUNX2','UGCG','PDE4B','PTPRS','BCL11A','PLXNA4'],
#                'Kaminsky@cDC1':['CLNK','CPVL','ASAP1','CCSER1','CADM1','TOX','CCDC26','CAMK2D','NCALD'],
#                'Kaminsky@DC_Mature':['REL','BIRC3','ARAP2','FAM49A','ARHGAP10','RNF145','CCR7','CCL22','CD86'],
#                'Kaminsky@DC_Langerhans':['DPYD','IRAK3','ARHGAP26','ITPR2','ETS2','PRKCE']}
# sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)
# sctri.plot_heterogeneity('Sun','dendritic_cells','violin',subset=subset,genes=['CADM1'])

# sc.pl.umap(sctri.adata,color=['ITGAM','ITGAX','PLTP','LILRB5','FCGR1A','MERTK','ADGRE1','HLA-DRA','HLA-DRB5','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DQA2','HLA-DQB2','HLA-DPA1','HLA-DPB1','CCR2','IGSF21'],cmap=bg_greyed_cmap('viridis'),vmin=1e-5)
# plt.savefig('output_opt1_RN_two_latest/reviewer_genes.pdf',bbox_inches='tight');plt.close()

# highlight4: monocyte
subset = ['cMonocyte','ncMonocyte']
sctri.plot_heterogeneity('Sun','monocytes','heatmap',col='Kaminsky',subset=subset)
for gene in ['CD14','FCGR3A']:
    sctri.plot_heterogeneity('Sun','monocytes','single_gene',col='Kaminsky',subset=subset,single_gene=gene)
sctri.plot_heterogeneity('Sun','monocytes','dual_gene',col='Kaminsky',subset=subset,dual_gene=['CD14','FCGR3A'])
    
























