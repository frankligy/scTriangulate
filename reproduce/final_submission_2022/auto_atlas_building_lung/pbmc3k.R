library(Seurat)
library(dplyr)
library(patchwork)

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")  # dgCMatrix S4
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) # Seurat S4

"assays, list of Assay object (from SeuratObject package)
 meta.data, data.frame
 active.assay, str vector
 active.ident, a factor indicating which column in meta.data to use as the cluster labels
 graph, list
 neighbors, list
 reductions, list
 images, list
 project.name, str vector
 version, S3 object, or list without name
 commands, list
 tools list
"

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # have column percent.mt in pbmc@meta.data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # pbmc@assay[['RNA']]@data will change
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # pbmc@assay[['RNA']]@meta.features and @var.features will change
all.genes <- rownames(pbmc);pbmc <- ScaleData(pbmc, features = all.genes) # pbmc@assay[['RNA']]@scale.data will change
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) # pbmc@reductions list will add $pca, which is a DimReduc S4 object from SeuratObject package
pbmc <- FindNeighbors(pbmc, dims = 1:10) # pbmc@graph add $RNA_nn and $RNA_snn, both are Graph S4 object from SeuratObject pacakge
pbmc <- FindClusters(pbmc, resolution = 0.5)  # will add cluster to pbmc@meta.data, also active.ident will automatically change
pbmc <- RunUMAP(pbmc, dims = 1:10) # pbmc@reduction list will add @umap 
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)

# essential command in Seurat: https://satijalab.org/seurat/articles/essential_commands.html

