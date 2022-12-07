setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/multiome_pbmc_qc')

library(Seurat)
library(cowplot)
library(dplyr)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(SeuratDisk)


# new analysis starting with h5ad files after QC
Convert('./post_rna_qc.h5ad','post_rna_qc.h5seurat',assay='RNA')
Convert('./post_atac_qc.h5ad','post_atac_qc.h5seurat',assay='ATAC')

seurat_rna <- LoadH5Seurat('./post_rna_qc.h5seurat')
seurat_atac <- LoadH5Seurat('./post_atac_qc.h5seurat')

seurat_multi = seurat_rna
seurat_multi@assays$ATAC = seurat_atac@assays$ATAC

DefaultAssay(seurat_multi) <- 'RNA'
seurat_multi <- NormalizeData(seurat_multi) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()


DefaultAssay(seurat_multi) <- 'ATAC'
seurat_multi <- RunTFIDF(seurat_multi)
seurat_multi <- FindTopFeatures(seurat_multi, min.cutoff = 'q0')
seurat_multi <- RunSVD(seurat_multi)

seurat_multi <- FindMultiModalNeighbors(seurat_multi,reduction.list=list('pca','lsi'),
                        dims.list=list(1:50,2:50))
seurat_multi <- RunUMAP(seurat_multi,nn.name='weighted.nn',reduction.name='wnn.umap',reduction.key='wnnUMAP_')

for(r in seq(from=0.25,to=3.75,by=0.25)) {
  seurat_multi <- FindClusters(seurat_multi,graph.name='wsnn',algorithm=3,resolution=r)
}  


meta.data <- seurat_multi@meta.data
wnn.umap <- seurat_multi@reductions$wnn.umap@cell.embeddings

write.table(meta.data,file='wnn_metadata.txt',sep='\t',col.names = NA)
write.table(wnn.umap,file='wnn_umap.txt',sep='\t',col.names=NA)







