#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os
import sys
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap

'''
all the input file, intermediate results are deposited on https://www.synapse.org/#!Synapse:syn26320731
'''

# large_txt_to_mtx(int_file='GSE161382_counts_matrix.txt',out_folder='./data',gene_is_index=True,type_convert_to='int16') 
adata = mtx_to_adata(int_folder='./data',gene_is_index=True,feature='genes')
adata = just_log_norm(adata)
add_annotations(adata,'./Groups/groups.Columbia-round-2.txt',['label'],0,['Columbia'])
add_annotations(adata,'./Groups/groups.Kaminsky-round-2.txt',['label'],0,['Kaminsky'])
add_annotations(adata,'./Groups/groups.Krasnow-round-2.txt',['label'],0,['Krasnow'])
add_annotations(adata,'./Groups/groups.Rhesus-round-2.txt',['label'],0,['Rhesus'])
add_annotations(adata,'./Groups/groups.Sun.txt',['label'],0,['Sun'])
add_umap(adata,'GSE161382_UMAP_coord.txt','pandas',['umap_x','umap_y'],0)
adata.write('input.h5ad')
umap_dual_view_save(adata,cols=['Columbia','Kaminsky','Krasnow','Rhesus','Sun'])

adata = sc.read('input.h5ad')

# run
sctri = ScTriangulate(dir='output_opt1_RN_two_new',adata=adata,query=['Sun','Columbia','Kaminsky','Krasnow'])
sctri.lazy_run(win_fraction_cutoff=0.4)


# insights
sctri = ScTriangulate.deserialize('output_opt1_RN_two/after_pruned_assess.p')
sctri.adata.obs['pruned_final_annotation'] = [item.split('@')[0] for item in sctri.adata.obs['pruned']]
sctri.plot_umap('pruned_final_annotation','category')
for col in ['Sun','Krasnow','Kaminsky','Columbia','Rhesus','final_annotation','pruned']:    
    sctri.plot_umap(col,'category',format='png')
sctri.plot_umap('confidence','continuous',umap_cmap='viridis')


# hightlight1: AT1
subset = ['Columbia@AT1','Kaminsky@ATI']
sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','umap',subset=subset,format='png')
sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','violin',subset=subset,format='pdf',genes=['IGFBP2'])
marker_dict = sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','heatmap',subset=subset)
marker_dict['Kaminsky@ATI'].to_csv('output_opt1_RN_two/Kaminsky@ATI.txt',sep='\t')
sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','violin',subset=subset,format='pdf',genes=['IGF2BP2'])
sctri.plot_heterogeneity('Sun','alveolar_type_1_cells','single_gene',subset=subset,format='pdf',single_gene='IGF2BP2')

# hightlight2: Alveolar_macrophage
subset = ['Kaminsky@Macrophage_Alveolar','Sun@Alveolar_macrophages']
marker_dict = sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap',subset=subset,format='pdf')
marker_dict['Kaminsky@Macrophage_Alveolar'].to_csv('output_opt1_RN_two/Kaminsky.txt',sep='\t')
marker_dict['Sun@Alveolar_macrophages'].to_csv('output_opt1_RN_two/Sun.txt',sep='\t')
sctri.plot_heterogeneity('Sun','Alveolar_macrophages','umap',subset=subset,format='png')
sctri.plot_heterogeneity('Sun','Alveolar_macrophages','violin',subset=subset,format='pdf',genes=['INHBA','FABP4'])
sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='INHBA')
sctri.plot_heterogeneity('Sun','Alveolar_macrophages','single_gene',subset=subset,format='pdf',single_gene='FABP4')
marker_gene_dict = {'Kaminsky@Macrophage_Alveolar':['C1QB','C1QA','C1QC','INHBA','FABP4','HLA-DPA1','CD52','HLA-DPA1','HLA-DRB1'],
                    'Sun@Alveolar_macrophages':['DOCK4','HIF1A','VCAN','ADAM9','STARD13','PLPP3','MERTK','GNA12']}
sctri.plot_heterogeneity('Sun','Alveolar_macrophages','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)

# highlight3: dendritic cells
subset = ['Kaminsky@cDC1','Krasnow@Plasmacytoid_Dendritic','Kaminsky@DC_Mature','Kaminsky@DC_Langerhans']
marker_dict = sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap',subset=subset,format='pdf')
marker_dict['Kaminsky@cDC1'].to_csv('output_opt1_RN_two/Kaminsky@cDC1.txt',sep='\t')
marker_dict['Kaminsky@DC_Mature'].to_csv('output_opt1_RN_two/Kaminsky@DC_Mature',sep='\t')
marker_dict['Kaminsky@DC_Langerhans'].to_csv('output_opt1_RN_two/Kaminsky@DC_Langerhans',sep='\t')
marker_dict['Krasnow@Plasmacytoid_Dendritic'].to_csv('output_opt1_RN_two/Krasnow@Plasmacytoid_Dendritic',sep='\t')
sctri.plot_heterogeneity('Sun','dendritic_cells','umap',format='png',subset=subset)
sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap',subset=subset)
marker_gene_dict = {'Krasnow@Plasmacytoid_Dendritic':['TCF4','ZFAT','CCDC50','RUNX2','UGCG','PDE4B','PTPRS','BCL11A','PLXNA4'],
               'Kaminsky@cDC1':['CLNK','CPVL','ASAP1','CCSER1','CADM1','TOX','CCDC26','CAMK2D','NCALD'],
               'Kaminsky@DC_Mature':['REL','BIRC3','ARAP2','FAM49A','ARHGAP10','RNF145','CCR7','CCL22','CD86'],
               'Kaminsky@DC_Langerhans':['DPYD','IRAK3','ARHGAP26','ITPR2','ETS2','PRKCE']}
sctri.plot_heterogeneity('Sun','dendritic_cells','heatmap+umap',subset=subset,marker_gene_dict=marker_gene_dict)
sctri.plot_heterogeneity('Sun','dendritic_cells','violin',subset=subset,genes=['CADM1'])






















