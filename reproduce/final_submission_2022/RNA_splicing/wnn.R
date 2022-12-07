setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/CPX-scRNA-Seq_Splicing')

library(Seurat)
library(cowplot)
library(dplyr)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(SeuratDisk)

Convert('./to_wnn_rna.h5ad','./wnn_rna.h5seurat',assay='RNA')
Convert('./to_wnn_splice.h5ad','./wnn_splice.h5seurat',assay='SPLICE')

seurat_rna = LoadH5Seurat('./wnn_rna.h5seurat')
seurat_splice = LoadH5Seurat('./wnn_splice.h5seurat')

seurat_combine = seurat_rna
seurat_combine@assays$SPLICE = seurat_splice@assays$SPLICE

DefaultAssay(seurat_combine) <- 'RNA'
seurat_combine <- NormalizeData(seurat_combine) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(seurat_combine) <- 'SPLICE'
seurat_combine <- FindVariableFeatures(seurat_combine) %>% ScaleData() %>% RunPCA(reduction.name = 'spca')

seurat_combine <- FindMultiModalNeighbors(seurat_combine,reduction.list=list('pca','spca'),
                                      dims.list=list(1:50,1:30))
seurat_combine <- RunUMAP(seurat_combine,nn.name='weighted.nn',reduction.name='wnn.umap',reduction.key='wnnUMAP_')

for(r in seq(from=0.25,to=3.75,by=0.25)) {
  seurat_combine <- FindClusters(seurat_combine,graph.name='wsnn',algorithm=3,resolution=r)
}  

Idents(seurat_combine) = seurat_combine$wsnn_res.1.5
r1.5.splice.markers = FindAllMarkers(seurat_combine,assay='SPLICE')
write.table(r1.5.splice.markers,file='wnn_splice_marker_r1.5.txt',sep='\t',col.names=NA)
r1.5.gene.markers = FindAllMarkers(seurat_combine,assay='RNA')
write.table(r1.5.gene.markers,file='wnn_gene_marker_r1.5.txt',sep='\t',col.names=NA)



meta.data <- seurat_combine@meta.data
wnn.umap <- seurat_combine@reductions$wnn.umap@cell.embeddings

write.table(meta.data,file='wnn_metadata.txt',sep='\t',col.names = NA)
write.table(wnn.umap,file='wnn_umap.txt',sep='\t',col.names=NA)

p1 <- DimPlot(seurat_combine, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()




