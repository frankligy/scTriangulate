setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/TNC1_qc')

library(Seurat)
library(cowplot)
library(dplyr)



tnc2.h5 <- Read10X_h5('28WM_ND19-341__TNC-RNA-ADT.h5')
tnc2.gene.seurat <- CreateSeuratObject(counts = tnc2.h5$`Gene Expression`)
tnc2.adt.assay <- CreateAssayObject(counts = tnc2.h5$`Antibody Capture`)
tnc2.gene.seurat[['ADT']] <- tnc2.adt.assay

tnc2.gene.seurat <- NormalizeData(tnc2.gene.seurat) %>% FindVariableFeatures() %>% 
                  ScaleData() %>% RunPCA()
DefaultAssay(tnc2.gene.seurat) <- 'ADT'
VariableFeatures(tnc2.gene.seurat) <- rownames(tnc2.gene.seurat[['ADT']])
tnc2.gene.seurat <- NormalizeData(tnc2.gene.seurat,normalizetion.method = 'CLR',
                                  margin = 2) %>% ScaleData() %>% 
                       RunPCA(reduction.name = 'apca')

tnc2.gene.seurat <- FindMultiModalNeighbors(tnc2.gene.seurat,
                                            reduction.list = list('pca','apca'),
                                            dims.list = list(1:30,1:18),
                                            modality.weight.name = 'RNA.weight')

tnc2.gene.seurat <- RunUMAP(tnc2.gene.seurat,nn.name='weighted.nn',
                            reduction.name='wnn.umap',reduction.key='wnnuMAP_')
for(r in seq(from=0.25,to=3.75,by=0.25)) {
tnc2.gene.seurat <- FindClusters(tnc2.gene.seurat,graph.name='wsnn',algorithm = 3,
                                 resolution=r, verbose=FALSE)
}
# p1 <- DimPlot(tnc2.gene.seurat,reduction='wnn.umap',label=TRUE,repel=TRUE,
#               label.size=2.5) + NoLegend()


meta.data <- tnc2.gene.seurat@meta.data
wnn.umap <- tnc2.gene.seurat@reductions$wnn.umap@cell.embeddings

write.table(meta.data,file='wnn_metadata.txt',sep='\t',col.names = NA)
write.table(wnn.umap,file='wnn_umap.txt',sep='\t',col.names=NA)








