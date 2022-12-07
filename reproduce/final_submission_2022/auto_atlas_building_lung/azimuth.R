library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)


ref <- pbmc

colormap <- list(annotation = CreateColorMap(object = ref, seed = 2))
colormap[["annotation"]] <- colormap[["annotation"]][sort(x = names(x = colormap[["annotation"]]))]

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("seurat_clusters"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

ref.dir = "/Users/ligk2e/R_envs/seurat"
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))