setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/TEA_qc')

library(Seurat)
library(cowplot)
library(dplyr)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(SeuratDisk)


# new analysis starting with h5ad files after QC
Convert('./post_rna_qc.h5ad','post_rna_qc.h5seurat',assay='RNA')
Convert('./post_adt_qc.h5ad','post_adt_qc.h5seurat',assay='ADT')
Convert('./post_atac_qc.h5ad','post_atac_qc.h5seurat',assay='ATAC')

seurat_rna <- LoadH5Seurat('./post_rna_qc.h5seurat')
seurat_adt <- LoadH5Seurat('./post_adt_qc.h5seurat')
seurat_atac <- LoadH5Seurat('./post_atac_qc.h5seurat')

seurat_tea = seurat_rna
seurat_tea@assays$ATAC = seurat_atac@assays$ATAC
seurat_tea@assays$ADT = seurat_adt@assays$ADT

DefaultAssay(seurat_tea) <- 'RNA'
seurat_tea <- NormalizeData(seurat_tea) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(seurat_tea) <- 'ADT'
VariableFeatures(seurat_tea) <- rownames(seurat_tea[["ADT"]])
seurat_tea <- NormalizeData(seurat_tea, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

DefaultAssay(seurat_tea) <- 'ATAC'
seurat_tea <- RunTFIDF(seurat_tea)
seurat_tea <- FindTopFeatures(seurat_tea, min.cutoff = 'q0')
seurat_tea <- RunSVD(seurat_tea)

seurat_tea <- FindMultiModalNeighbors(seurat_tea,reduction.list=list('pca','apca','lsi'),
                        dims.list=list(1:50,1:18,2:50))
seurat_tea <- RunUMAP(seurat_tea,nn.name='weighted.nn',reduction.name='wnn.umap',reduction.key='wnnUMAP_')

for(r in seq(from=0.25,to=3.75,by=0.25)) {
  seurat_tea <- FindClusters(seurat_tea,graph.name='wsnn',algorithm=3,resolution=r)
}  


meta.data <- seurat_tea@meta.data
wnn.umap <- seurat_tea@reductions$wnn.umap@cell.embeddings

write.table(meta.data,file='wnn_metadata.txt',sep='\t',col.names = NA)
write.table(wnn.umap,file='wnn_umap.txt',sep='\t',col.names=NA)







