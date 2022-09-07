library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)

setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/lung_new_sun/Frank')
# ref <- readRDS('2021_MNP_Verse.RDS')
# ref.data <- ref@assays[['RNA']]@counts
# new_ref <- CreateSeuratObject(counts = ref.data, project = "MNP", min.cells = 3, min.features = 200)
# cells <- colnames(new_ref@assays[['RNA']]@counts)
# new_ref@meta.data <- ref@meta.data[cells,] # make sure the order is the same
# new_ref@meta.data$Clusters <- as.factor(new_ref@meta.data$Clusters)
# new_ref@active.ident <- new_ref@meta.data$Clusters

# new_ref <- PercentageFeatureSet(new_ref, pattern = "^MT-", col.name = "percent.mt")
# new_ref <- SCTransform(new_ref, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
# new_ref <- RunPCA(new_ref, verbose = TRUE)
# new_ref <- RunUMAP(new_ref, dims = 1:30, verbose = TRUE,umap.method='uwot',return.model = TRUE)
# new_ref <- FindNeighbors(new_ref, dims = 1:30, verbose = TRUE)

# saveRDS(new_ref,'tmp_new_ref.rds')

# ref <- readRDS('tmp_new_ref.rds')

# colormap <- list(Clusters = CreateColorMap(object = ref, seed = 2))
# colormap[["Clusters"]] <- colormap[["Clusters"]][sort(x = names(x = colormap[["Clusters"]]))]

# ref <- AzimuthReference(
#   object = ref,
#   refUMAP = "umap",
#   refDR = "pca",
#   refAssay = "SCT",
#   metadata = c("Clusters"),
#   dims = 1:50,
#   k.param = 31,
#   colormap = colormap,
#   reference.version = "1.0.0"
# )

# ref.dir = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/lung_new_sun/Frank"
# SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
# saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))


# run azimuth
query <- RunAzimuth(query='input_raw_count.h5ad',
                    reference='./',
                    annotation.levels='Clusters')

saveRDS(query,'query.rds')

# save result
query <- readRDS('query.rds')
write.table(query@meta.data,'MNP_mapping_results.txt',sep='\t',col.names=NA)



