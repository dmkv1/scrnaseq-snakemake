threads <- snakemake@params[["threads"]]
random_seed <- snakemake@params[["random_seed"]]
sample_name <- snakemake@params[["sample"]]
knn <- snakemake@params[["KNN"]]

input_sce <- snakemake@input[["sce"]]
output_sce <- snakemake@output[["sce_annotated"]]
output_plot_umaps <- snakemake@output[["plot_UMAPs"]]
output_markers <- snakemake@output[["markers"]]

# Logging
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Parallelization
library(BiocParallel)
bpp <- MulticoreParam(workers = threads)

# Load the SCE
sce <- readRDS(input_sce)

library(ggpubr)
library(scuttle)
library(scater)
library(scran)

# Normalization
set.seed(random_seed)
clusters_igraph <- quickCluster(sce,
                                method = "igraph",
                                BPPARAM = bpp)
sce$clusters_igraph <- clusters_igraph

set.seed(random_seed)
sce <- computeSumFactors(sce,
                         cluster = clusters_igraph,
                         min.mean = 0.1,
                         BPPARAM = bpp)
sce <- logNormCounts(sce,
                     size.factors = sce$sizeFactor,
                     BPPARAM = bpp)

# PCA
variance_decomposition <- modelGeneVar(sce, BPPARAM = bpp)
top_HVG <- getTopHVGs(variance_decomposition, fdr.threshold = 0.25) # nolint
set.seed(random_seed)
sce <- denoisePCA(sce, subset.row = top_HVG, technical = variance_decomposition)
percent_var <- attr(reducedDim(sce), "percentVar")

# chose a number of PCA dimensions which explain more than 0.5% variance
chosen_pc_num <-
  percent_var[which(percent_var > 0.5)] %>% length()
print(
  paste(
    "Calculated PCA, number of dimensions explaining more than 0.5% variance:",
    chosen_pc_num
  )
)

# UMAP
set.seed(random_seed)
sce <- runUMAP(
  sce,
  pca = chosen_pc_num,
  n_neighbors = 30,
  min_dist = 0.2,
  dimred = "PCA",
  BPPARAM = bpp
)

# Clustering
set.seed(random_seed)
nn_clust <- clusterCells(
  sce,
  use.dimred = "PCA",
  full = TRUE,
  BLUSPARAM = bluster::SNNGraphParam(
    k = knn,
    type = "rank",
    cluster.fun = "louvain",
    BPPARAM = bpp
  )
)
sce[[paste0("clusters_louvain_K", knn)]] <- nn_clust$clusters

# Immunological cell type annotation
library(celldex)
ref <- DatabaseImmuneCellExpressionData()

rownames(sce) <- rowData(sce)[["Symbol"]]

library(SingleR)
set.seed(random_seed)
predicted_types_main <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  de.method = "classic",
  BPPARAM = bpp
)
sce$cell_type_main <- predicted_types_main$labels

set.seed(random_seed)
predicted_types_fine <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.fine,
  de.method = "classic",
  BPPARAM = bpp
)
sce$cell_type_fine <- predicted_types_fine$labels

plot_umaps <- ggarrange(
  plotUMAP(
    sce,
    color_by = paste0("clusters_louvain_K", knn),
    text_by = paste0("clusters_louvain_K", knn)
  ),
  plotUMAP(
    sce,
    color_by = "cell_type_main",
    text_by = paste0("clusters_louvain_K", knn)
  ),
  ncol = 2, nrow = 1
)

# Cluster markers
library(Seurat)

rownames(sce) <- rowData(sce)[["ID"]]
seu <- as.Seurat(sce)
seu <- RenameAssays(seu, originalexp = "RNA")

Idents(seu) <- paste0("clusters_louvain_K", knn)
markers <- FindAllMarkers(
  seu,
  test.use = "wilcox",
  min.pct = 0.1, # expressed in at least 10% of the cells
  logfc.threshold = 0.25, # moderate logFC threshold
  min.diff.pct = 0.1, # difference in cells expressing gene between clusters
  assay = "RNA",
  slot = "data"
)

markers$symbol <- rowData(sce)[["Symbol"]][match(markers$gene, rowData(sce)[["ID"]])] # nolint

# Export outputs
saveRDS(sce, file = output_sce)
saveRDS(markers, file = output_markers)

## Export plots
ggsave(
  filename = output_plot_umaps,
  plot = plot_umaps,
  bg = "white",
  dpi = 300,
  units = "px",
  width = 2000,
  height = 1000
)

# Close log connection
sink(type = "message")
sink(type = "output")
close(log)
