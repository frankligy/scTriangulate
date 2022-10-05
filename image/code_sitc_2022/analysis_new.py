#!/opt/anaconda3/envs/sctri_env/bin/python3.8

import pandas as pd
import squidpy as sq
import scanpy as sc
import os,sys
from sctriangulate.colors import *
from sctriangulate.preprocessing import *
from sctriangulate import *
from sctriangulate.spatial import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import f1_score, confusion_matrix, ConfusionMatrixDisplay
from tqdm import tqdm
from functools import reduce
from scipy.sparse import csr_matrix
os.chdir('/Users/ligk2e/Desktop/sitc')
sys.path.insert(0,'')
from utils import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

sctriangulate_setting(backend='Agg')


'''
scRNA dataset
'''

# make sense of the scRNA data
# adata_sc = mtx_to_adata(int_folder='Wu_etal_2021_BRCA_scRNASeq',gene_is_index=True,feature='count_matrix_genes.tsv',feature_col='index',
#                         barcode='count_matrix_barcodes.tsv',barcode_col='index',matrix='count_matrix_sparse.mtx')
# add_annotations(adata_sc,inputs='Wu_etal_2021_BRCA_scRNASeq/metadata.csv',cols_input=['orig.ident','subtype','celltype_subset','celltype_minor','celltype_major'],sep=',')
# adata_sc.X = csr_matrix(adata_sc.X)
# adata_sc = make_sure_adata_writable(adata_sc)
# adata_sc.write('adata_sc.h5ad')

adata_sc = sc.read('adata_sc.h5ad')
# adata_sc_processed = adata_sc.copy()
# adata_sc_processed = scanpy_recipe(adata_sc_processed,species='human',is_log=False,resolutions=[1],modality='rna',pca_n_comps=50,n_top_genes=3000)
adata_sc_processed = sc.read('adata_after_scanpy_recipe_rna_1_umap_True.h5ad')
# umap_dual_view_save(adata_sc_processed,cols=['orig.ident','subtype','celltype_subset','celltype_minor','celltype_major'])

final_annotation = []
for major,minor in zip(adata_sc_processed.obs['celltype_major'],adata_sc_processed.obs['celltype_minor']):
    if minor == 'Macrophage' or minor == 'DCs' or minor == 'NK cells':
        final_annotation.append(minor)
    elif major == 'B-cells' or major == 'Plasmablasts':
        final_annotation.append('B cells')
    elif major == 'CAFs':
        final_annotation.append('fibroblast')
    elif major == 'Cancer Epithelial':
        final_annotation.append('Cancer Epithelial')
    elif major == 'Endothelial':
        final_annotation.append('Endothelial')
    elif major == 'Normal Epithelial':
        final_annotation.append('Normal Epithelial')
    elif major == 'T-cells':
        final_annotation.append('T cells')
    else:
        final_annotation.append('unknown')
adata_sc_processed.obs['final_annotation'] = final_annotation
adata_sc_processed = adata_sc_processed[adata_sc_processed.obs['final_annotation']!='unknown',:] # 91850 * 29733
adata_sc = adata_sc[adata_sc_processed.obs_names,:]  # 91850 * 29733
add_annotations(adata_sc,inputs=adata_sc_processed.obs,cols_input=['final_annotation'],kind='memory')
# umap_dual_view_save(adata_sc_processed,cols=['final_annotation'])

adata_sc_processed.uns['log1p']["base"] = None
sc.tl.rank_genes_groups(adata_sc_processed,'final_annotation')  # plot both annotation and coarse annotation
sc.tl.filter_rank_genes_groups(adata_sc_processed)
sc.pl.rank_genes_groups_heatmap(adata_sc_processed, key='rank_genes_groups_filtered',swap_axes=True)
# plt.savefig('final_annotation_marker_genes.pdf',bbox_inches='tight');plt.close()
markers = pd.DataFrame.from_records(adata_sc_processed.uns['rank_genes_groups_filtered']['names'])


'''
marker dict
'''

markers_dict = {
    'B cells':['CD79A','RALGPS2','CD79B','MS4A1','BANK1','TNFRSF13C','IGHM','MEF2C','CD19','CD22','CR2','FCRL2','PAX5'],
    'DCs':['MS4A6A','GPR183','HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DQB1','HLA-DQA1','HLA-DMA','CST3','IRF7'],
    'Normal Epithelial':['KRT7','PIGR','ELF3','CYB5A','KRT8','KRT19','TACSTD2','EPCAM'],
    'Cancer Epithelial':['KRT7','PIGR','ELF3','CYB5A','KRT8','KRT19','TACSTD2','EPCAM'],
    'fibroblast': ['BGN', 'SPARC', 'CALD1', 'COL6A1', 'COL1A2', 'C1S'],
    'Macrophage':['C1QB','C1QA','CTSZ','CTSL','FCER1G','C1QC','LYZ','CST1','CD68','FCN1','VCAN','MARCO','APOE'],
    'NK cells': ['NKG7', 'KLRD1','GNLY','KLRF1','NCR1','PTGDR','SH2D1B'],
    'T cells': ['CD3E','TRBC2','TRAC','CD2', 'CD3D','CD28','CD3G','CD5','CD6','PBX4','THEMIS','TRAT1'],
    'Endothelial':['CLDN5','AQP1','PECAM1','NPDC1','VWF','GNG11','RAMP2','IGFBP7','CLEC14A']
}


'''
spatial datasets
'''

# # make sense of the spatial data
# ids = ['1160920F','1142243F']  # ['CID44971','CID4535','CID4465','CID4290','1160920F','1142243F']
# for id_ in ids:
#     adata_spatial = read_spatial_data(mtx_folder='filtered_count_matrices/{}_filtered_count_matrix'.format(id_),
#                                       spatial_folder='filtered_count_matrices/{}_filtered_count_matrix/{}_spatial'.format(id_,id_),
#                                       spatial_library_id=id_,feature='features.tsv')
#     sc.pp.calculate_qc_metrics(adata_spatial, percent_top=None, inplace=True, log1p=False)
#     sc.pl.spatial(adata_spatial,library_id=id_,title=id_,alpha_img=1)
#     print(adata_spatial.shape[0])


'''
CID44971 matched
'''
adata_sc_processed_subset = adata_sc_processed[adata_sc_processed.obs['orig.ident']=='CID44971',:].copy()
adata_sc_subset = adata_sc[adata_sc.obs['orig.ident']=='CID44971',:].copy()
# umap_dual_view_save(adata_sc_processed_subset,cols=['final_annotation'])
id_ = 'CID44971'
adata_spatial = read_spatial_data(mode_count='mtx',mode_spatial='visium',
                                  mtx_folder='spatial/{}_filtered_count_matrix'.format(id_),
                                  spatial_folder='spatial/{}_filtered_count_matrix/{}_spatial'.format(id_,id_), spatial_library_id=id_,
                                  feature='features.tsv')
sc.pp.calculate_qc_metrics(adata_spatial, percent_top=None, inplace=True, log1p=False)
sc.pl.spatial(adata_spatial,library_id=id_,title=id_,alpha_img=1)
# plt.savefig('CID44971_spatial_image.pdf',bbox_inches='tight')
# plt.close()

# # prepare inputs and run
# prepare_inputs(adata_sc_subset,adata_sc_processed_subset,adata_spatial,name='CID44971',column='final_annotation')
# run_ingest(adata_sc_subset,adata_spatial,'final_annotation','inputs/decon_results')
# run_scanorama(adata_sc_subset,adata_spatial,'final_annotation','inputs/decon_results')


'''
CID4535 to CID44971
'''
# adata_sc_processed_subset = adata_sc_processed[adata_sc_processed.obs['orig.ident']=='CID4535',:].copy()
# adata_sc_subset = adata_sc[adata_sc.obs['orig.ident']=='CID4535',:].copy()
# umap_dual_view_save(adata_sc_processed_subset,cols=['final_annotation'])
# id_ = 'CID44971'
# adata_spatial = read_spatial_data(mtx_folder='spatial/{}_filtered_count_matrix'.format(id_),
#                                   spatial_folder='spatial/{}_filtered_count_matrix/{}_spatial'.format(id_,id_),
#                                   spatial_library_id=id_,feature='features.tsv')

# prepare_inputs(adata_sc_subset,adata_sc_processed_subset,adata_spatial,name='CID4535',column='final_annotation')
# run_ingest(adata_sc_subset,adata_spatial,'final_annotation','inputs/decon_results')
# run_scanorama(adata_sc_subset,adata_spatial,'final_annotation','inputs/decon_results')


'''
post-processing
'''

# ### deconRNASeq
# outdir = 'inputs/decon_results'
# decon_rna_prop = pd.read_csv(os.path.join(outdir,'DeconRNASeq_prop.txt'),sep='\t',index_col=0)
# decon_rna_prop.columns = [item.replace('.',' ') for item in decon_rna_prop.columns]
# decon_rna_prop.index = [item.replace('.','-') for item in decon_rna_prop.index]
# decon_rna_prop.to_csv(os.path.join(outdir,'DeconRNASeq_prop_processed.txt'),sep='\t')

# ### cell2location
# decon_rna_prop = pd.read_csv(os.path.join(outdir,'cell2location_prop.txt'),sep='\t',index_col=0)
# decon_rna_prop.columns = [item.replace('q05cell_abundance_w_sf_','') for item in decon_rna_prop.columns]
# decon_rna_prop = decon_rna_prop.apply(func=lambda x:x/x.sum(),axis=1,result_type='expand')
# decon_rna_prop.to_csv(os.path.join(outdir,'cell2location_prop_processed.txt'),sep='\t')

'''
evaluation
'''

## generate some plots for illustrating how the benchmark performance was evaluated
# sc.pl.spatial(adata_spatial,library_id=id_,alpha_img=1,color=['C1QB','C1QA','CTSZ','CTSL','FCER1G','C1QC','LYZ','CST1','CD68','FCN1','VCAN','MARCO','APOE'])
# plt.savefig('CID44971_macrophage.pdf',bbox_inches='tight')
# plt.close()

# sc.pl.spatial(adata_spatial,library_id=id_,alpha_img=1,color=['CD3E','TRBC2','TRAC','CD2', 'CD3D','CD28','CD3G','CD5','CD6','PBX4','THEMIS','TRAT1'])
# plt.savefig('CID44971_Tcell.pdf',bbox_inches='tight')
# plt.close()

# for gene in ['CD3E','CD3D','CD2']:
#     evaluate_decon_result_single('inputs/decon_results/cell2location_prop_processed.txt',adata_spatial,output=None,outdir='inputs/decon_results',plot_info=('T cells',gene))


# visualize_decon_result('inputs/decon_results/CARD_prop.txt',adata_spatial,output='CARD',outdir='inputs/decon_results',markers_to_plot=[])
# evaluate_decon_result_merge('inputs/decon_results/CARD_prop.txt',adata_spatial,output='CARD',outdir='inputs/decon_results',markers_dict=markers_dict)
# visualize_decon_result('inputs/decon_results/cell2location_prop_processed.txt',adata_spatial,output='cell2location',outdir='inputs/decon_results',markers_to_plot=[])
# evaluate_decon_result_merge('inputs/decon_results/cell2location_prop_processed.txt',adata_spatial,output='cell2location',outdir='inputs/decon_results',markers_dict=markers_dict)

# decon = pd.read_csv('inputs/decon_results/cell2location_prop_processed.txt',index_col=0,sep='\t')
# plot_deconvolution(adata_spatial,decon,fraction=1,alpha=0.3)

# tmp_df_list = []
# tmp_df_list.append(build_the_df_merge('ingest_un','inputs/decon_results/ingest_prop.txt',adata_spatial,0,markers_dict))
# tmp_df_list.append(build_the_df_merge('scanoroma_un','inputs/decon_results/scanoroma_prop.txt',adata_spatial,0,markers_dict))
# tmp_df_list.append(build_the_df_merge('seurat_un','inputs/decon_results/seurat_prop.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('DeconRNASeq_un','inputs/decon_results/DeconRNASeq_prop_processed.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('CARD_un','inputs/decon_results/CARD_prop.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('cell2location_un','inputs/decon_results/cell2location_prop_processed.txt',adata_spatial,0.05,markers_dict))

# tmp_df_list.append(build_the_df_merge('un_ingest','decon_results_A/ingest_prop.txt',adata_spatial,0,markers_dict))
# tmp_df_list.append(build_the_df_merge('un_scanoroma','decon_results_A/scanoroma_prop.txt',adata_spatial,0,markers_dict))
# tmp_df_list.append(build_the_df_merge('un_seurat','decon_results_A/seurat_prop.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('un_DeconRNASeq','decon_results_A/DeconRNASeq_prop_processed.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('un_CARD','decon_results_A/CARD_prop.txt',adata_spatial,0.05,markers_dict))
# tmp_df_list.append(build_the_df_merge('un_cell2location','decon_results_A/cell2location_prop_processed.txt',adata_spatial,0.05,markers_dict))
# metrics = pd.concat(tmp_df_list,axis=0).fillna(0)

# for metric in ['f1_score','spearmanr','pearsonr']:
#     fig,ax = plt.subplots()
#     sns.barplot(x='cell_type',y=metric,data=metrics,hue='method',palette='tab20',ax=ax)
#     plt.savefig('{}_per_cell_type.pdf'.format(metric),bbox_inches='tight')
#     plt.close()

# for metric in ['f1_score','spearmanr','pearsonr']:
#     # combine cell types
#     data = []
#     score = metric
#     for method,sub_df in metrics.groupby(by='method'):
#         data.append((method,sub_df[score].mean()))
#     metrics_combine = pd.DataFrame.from_records(data,columns=['method',score]).sort_values(by=score)
#     fig,ax = plt.subplots()
#     sns.barplot(x='method',y=score,data=metrics_combine,palette='tab20',ax=ax)
#     plt.savefig('{}_whole.pdf'.format(metric),bbox_inches='tight')
#     plt.close()



# spatial structure
id_ = 'CID44971'
adata_spatial = read_spatial_data(mtx_folder='spatial/{}_filtered_count_matrix'.format(id_),
                                  spatial_folder='spatial/{}_filtered_count_matrix/{}_spatial'.format(id_,id_),
                                  spatial_library_id=id_,feature='features.tsv')
decon = pd.read_csv('inputs/decon_results/cell2location_prop_processed.txt',sep='\t',index_col=0)
adata_spatial.obsm['decon'] = decon.loc[adata_spatial.obs_names,:]
# adata_neigh = identify_ecosystem(adata_spatial,n_rings=2)
# adata_neigh.write('adata_neigh.h5ad')

adata_neigh = sc.read('adata_neigh.h5ad')
for c in range(19):
    sc.pl.spatial(adata_neigh,color='leiden',groups=[str(c)],alpha_img=0.2)
    plt.savefig('{}_ecosystem_plot.pdf'.format(c),bbox_inches='tight')
    plt.close()
sys.exit('stop')
df_subgraph = identify_spatial_program(adata_spatial,n_rings=1)

# using sctriangulate
spatial_adata_coord = create_spatial_features(adata_spatial,mode='coordinate')
sc.pp.neighbors(spatial_adata_coord)
sc.tl.leiden(spatial_adata_coord)
spatial_adata_coord.uns['spatial'] = adata_spatial.uns['spatial']
spatial_adata_coord.obsm['spatial'] = adata_spatial.obsm['spatial']
sc.pl.spatial(spatial_adata_coord,color='leiden')

spatial_adata_graph = create_spatial_features(adata_spatial,mode='graph_importance')
sc.pp.neighbors(spatial_adata_graph)
sc.tl.leiden(spatial_adata_graph)
spatial_adata_graph.uns['spatial'] = adata_spatial.uns['spatial']
spatial_adata_graph.obsm['spatial'] = adata_spatial.obsm['spatial']
sc.pl.spatial(spatial_adata_graph,color='leiden')

spatial_adata_image = create_spatial_features(adata_spatial,mode='tissue_image',library_id='CID44971',img_key='hires',
                                              segmentation_feature=True)
sc.pp.neighbors(spatial_adata_image)
sc.tl.leiden(spatial_adata_image)
spatial_adata_image.uns['spatial'] = adata_spatial.uns['spatial']
spatial_adata_image.obsm['spatial'] = adata_spatial.obsm['spatial']
sc.pl.spatial(spatial_adata_image,color='leiden')

sc.pp.normalize_total(adata_spatial, target_sum=1e4)
sc.pp.log1p(adata_spatial)
sc.pp.highly_variable_genes(adata_spatial, flavor='seurat', n_top_genes=3000)
adata_spatial = adata_spatial[:, adata_spatial.var['highly_variable']]
sc.pp.pca(adata_spatial)
sc.pp.neighbors(adata_spatial)
sc.tl.leiden(adata_spatial)
sc.pl.spatial(adata_spatial,color='leiden')

### triangulate transcriptome and spatial_adata_coord
adata_spatial.obs.rename(columns={'leiden':'leiden_rna'},inplace=True)
spatial_adata_coord.obs.rename(columns={'leiden':'leiden_coord'},inplace=True)
adata_merge = concat_rna_and_other(adata_spatial,spatial_adata_coord,umap='rna',umap_key='spatial',name='spatial',prefix='spatial_')
sctri = ScTriangulate(dir='output',adata=adata_merge,query=['leiden_rna','leiden_coord'])
sctri.lazy_run(scale_sccaf=False,assess_pruned=True,viewer_cluster=True,viewer_heterogeneity=True)