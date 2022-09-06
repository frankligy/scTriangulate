library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)

setwd('/data/salomonis2/LabFiles/Frank-Li/scTriangulate/lung_new_sun/Frank')
ref <- readRDS('2021_MNP_Verse.RDS')
ref.data <- ref@assays[['RNA']]@counts
new_ref <- CreateSeuratObject(counts = ref.data, project = "MNP", min.cells = 3, min.features = 200)
cells <- colnames(new_ref@assays[['RNA']]@counts)
new_ref@meta.data <- ref@meta.data[cells,] # make sure the order is the same
new_ref@meta.data$Clusters <- as.factor(new_ref@meta.data$Clusters)
new_ref@active.ident <- new_ref@meta.data$Clusters

new_ref <- PercentageFeatureSet(new_ref, pattern = "^MT-", col.name = "percent.mt")
new_ref <- SCTransform(new_ref, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
new_ref <- RunPCA(new_ref, verbose = FALSE)
new_ref <- RunUMAP(new_ref, dims = 1:30, verbose = FALSE,umap.method='uwot',return.model = TRUE)
new_ref <- FindNeighbors(new_ref, dims = 1:30, verbose = FALSE)

saveRDS(new_ref,'tmp_new_ref.rds')
quit()

colormap <- list(annotation = CreateColorMap(object = ref, seed = 2))
colormap[["Cluster"]] <- colormap[["Cluster"]][sort(x = names(x = colormap[["Cluster"]]))]

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("Cluster"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

ref.dir = "/data/salomonis2/LabFiles/Frank-Li/scTriangulate/lung_new_sun/Frank"
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(ref.dir, "ref.Rds"))


# run Azimuth
# Load the reference
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "/reference-data/human_lung_v2/")

# Load the query object for mapping
# Change the file path based on where the query file is located on your system.
query <- LoadFileInput(path = "input.h5ad")
query <- ConvertGeneNames(
  object = query,
  reference.names = rownames(x = reference$map),
  homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
)

# Calculate nCount_RNA and nFeature_RNA if the query does not
# contain them already
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
  calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
  colnames(x = calcn) <- paste(
    colnames(x = calcn),
    "RNA",
    sep = '_'
  )
  query <- AddMetaData(
    object = query,
    metadata = calcn
  )
  rm(calcn)
}


# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

# Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
cells.use <- query[["nCount_RNA", drop = TRUE]] <= 4452 &
  query[["nCount_RNA", drop = TRUE]] >= 859 &
  query[["nFeature_RNA", drop = TRUE]] <= 3997 &
  query[["nFeature_RNA", drop = TRUE]] >= 501

# If the query contains mitochondrial genes, filter cells based on the
# thresholds for percent.mt you set in the app
if ("percent.mt" %in% c(colnames(x = query[[]]))) {
  cells.use <- cells.use & (query[["percent.mt", drop = TRUE]] <= 9 &
                              query[["percent.mt", drop = TRUE]] >= 0)
}

# Remove filtered cells from the query
query <- query[, cells.use]

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = "ann_finest_level", function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- "ann_finest_level"
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
id <- "ann_finest_level"[1]
predicted.id <- paste0("predicted.", id)

# DimPlot of the reference
DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()

# Plot the prediction score for the class CD16 Mono
FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
VlnPlot(object = query, features = "CD16 Mono", group.by = predicted.id) + NoLegend()

# Plot an RNA feature
FeaturePlot(object = query, features = "EMP2", reduction = "proj.umap")
VlnPlot(object = query, features = "EMP2", group.by = predicted.id) + NoLegend()

# Plot an imputed protein feature
if (FALSE) {
  FeaturePlot(object = query, features = "CD3-1", reduction = "proj.umap")
  VlnPlot(object = query, features = "CD3-1", group.by = predicted.id) + NoLegend()
}






