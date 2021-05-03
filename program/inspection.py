#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import scanpy as sc
import matplotlib.pyplot as plt

def plot_DE_umap_save(adata,reference):
    with open('./scTriangulate_inspection/log.txt','w') as f:
        for cluster in adata.obs[reference].astype('category').cat.categories:
            adata_s = adata[adata.obs[reference]==cluster,:].copy()

            # first, save adata_s.h5ad (cellxgene)
            adata_s.write('./scTriangulate_inspection/to_cellxgene_{}.h5ad'.format(cluster))

            # second, umap
            sc.pl.umap(adata_s,color=['reassign_prefix'])
            plt.savefig('./scTriangulate_inspection/UMAP_{}.pdf'.format(cluster),bbox_inches='tight')
            plt.close()

            # third, DE
            if adata_s.uns.get('rank_genes_groups') != None:
                del adata_s.uns['rank_genes_groups']
            if len(adata_s.obs['reassign_prefix'].unique()) == 1: # it is already unique
                print('{0} entirely being assigned to one type, no need to do DE'.format(cluster),file=f)
                continue
            sc.tl.rank_genes_groups(adata_s,groupby='reassign_prefix')
            number_of_groups = len(adata_s.obs['reassign_prefix'].unique())
            genes_to_pick = 50 // number_of_groups
            sc.pl.rank_genes_groups_heatmap(adata_s,n_genes=genes_to_pick,swap_axes=True)
            plt.savefig('./scTriangulate_inspection/DE_heatmap_{}.pdf'.format(cluster),bbox_inches='tight')
            plt.close()






        


