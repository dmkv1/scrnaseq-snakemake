threads <- snakemake@params[["threads"]]
random_seed <- snakemake@params[["random_seed"]]
sample_name <- snakemake@params[["sample"]]
mito_genes_perc_thresh <- snakemake@params[["mito_perc"]]
knn <- snakemake@params[["KNN"]]
cluster_removal_threshold <- snakemake@params[["cluster_perc"]]

input_sce <- snakemake@input[["sce"]]
output_sce_unfiltered <- snakemake@output[["sce_unfiltered"]]
output_sce_filtered <- snakemake@output[["sce_filtered"]]
output_qc_plot_transcripts <- snakemake@output[["plot_transcripts"]]
output_qc_plot_doublets <- snakemake@output[["plot_doublets"]]
output_qc_plot_umaps <- snakemake@output[["plot_UMAPs"]]

# Parallelization
library(BiocParallel)
bpp <- MulticoreParam(workers = threads)

# Logging
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load the SCE
sce <- readRDS(input_sce)

# Median absolute deviation (MAD) based adaptive QC
library(scuttle)
qc_genes <- list(Mito = stringr::str_which(rowData(sce)[["Symbol"]], "^MT-"))
sce <- addPerCellQCMetrics(sce, subsets = qc_genes)
qc_filter <- quickPerCellQC(colData(sce),
                            sub.fields = c("subsets_Mito_percent"))
colData(sce) <- cbind(colData(sce), qc_filter)

sce$subsets_Mito_filter_out <- sce$subsets_Mito_percent > mito_genes_perc_thresh

# QC vizualization
library(scater)
library(ggpubr)

plot_qc_transcripts <- ggarrange(
  plotColData(
    sce,
    x = "Sample",
    y = "sum",
    colour_by = "discard"
  ) +
    scale_y_log10() +
    labs(title = "Total transcripts per cell") +
    theme(panel.grid.major.y = element_line(
      colour = "gray70", linetype = "dashed"
    )),
  plotColData(
    sce,
    x = "Sample",
    y = "detected",
    colour_by = "discard"
  ) +
    scale_y_log10() +
    labs(title = "UMI per cell") +
    theme(panel.grid.major.y = element_line(
      colour = "gray70", linetype = "dashed"
    )),
  plotColData(
    sce,
    x = "Sample",
    y = "subsets_Mito_percent",
    colour_by = "discard"
  ) +
    geom_hline(
      yintercept = mito_genes_perc_thresh,
      color = "firebrick",
      linetype = "dashed"
    ) +
    labs(title = "% of mitochondrial genes") +
    theme(panel.grid.major.y = element_line(
      colour = "gray70", linetype = "dashed"
    )),
  plotColData(
    sce,
    x = "detected",
    y = "sum",
    colour_by = "discard"
  ) +
    scale_x_log10() +
    scale_y_log10() +
    theme(panel.grid.major = element_line(
      colour = "gray70", linetype = "dashed"
    )) +
    labs(y = "Total transcripts per cell", x = "UMI per cell"),
  ncol = 2,
  nrow = 2,
  common.legend = TRUE,
  legend = "right"
)

# Doublet detection
library(scDblFinder)
set.seed(random_seed)
sce <- scDblFinder(sce, clusters = FALSE)
plot_doublets <- plotColData(sce,
  x = "detected",
  y = "sum",
  colour_by = "scDblFinder.class"
) +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.grid.major =
          element_line(colour = "gray70", linetype = "dashed")) +
  labs(y = "Total transcripts per cell", x = "UMI per cell")

# Simplified normalization
library(scran)
sce <- logNormCounts(sce, BPPARAM = bpp)

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

plot_qc_umaps <- ggarrange(
  plotUMAP(
    sce,
    color_by = paste0("clusters_louvain_K", knn),
    text_by = paste0("clusters_louvain_K", knn)
  ),
  plotUMAP(
    sce,
    color_by = "subsets_Mito_filter_out",
    text_by = paste0("clusters_louvain_K", knn),
  ) +
    labs(color = "% mitochondrial genes"),
  plotUMAP(
    sce,
    color_by = "discard",
    text_by = paste0("clusters_louvain_K", knn)
  ),
  plotUMAP(
    sce,
    color_by = "scDblFinder.class",
    text_by = paste0("clusters_louvain_K", knn)
  ),
  ncol = 2,
  nrow = 2
)

# Remove clusters with more than <cluster_removal_threshold>% of "bad" cells
qc_data <- as.data.frame(colData(sce)[, c(
  paste0("clusters_louvain_K", knn),
  "discard",
  "subsets_Mito_filter_out"
)])
qc_data$both_filters <- qc_data$discard &
  qc_data$subsets_Mito_filter_out
colnames(qc_data) <- c(
  "cluster",
  "discard",
  "subsets_Mito_filter_out",
  "both_filters"
)

cluster_stats <- aggregate(
  both_filters ~ cluster,
  data = qc_data,
  FUN = function(x) {
    sum(x) / length(x)
  }
)
cluster_stats <-
  cluster_stats[order(cluster_stats$both_filters, decreasing = TRUE), ]

cat("Proportion of damaged cells per cluster:")
print(cluster_stats)

clusters_to_remove <-
  cluster_stats$cluster[cluster_stats$both_filters >
                        (cluster_removal_threshold / 100)]

original_count <- ncol(sce)

is_bad_qc <- colData(sce)$discard |
  colData(sce)$subsets_Mito_filter_out |
  colData(sce)$scDblFinder.class == "doublet"

if (length(clusters_to_remove) > 0) {
  is_bad_cluster <- colData(sce)[[paste0("clusters_louvain_K", knn)]] %in%
    clusters_to_remove
  keep_cells <- !(is_bad_cluster | is_bad_qc)

  sce_filtered <- sce[, keep_cells]

  removed_by_cluster <- sum(is_bad_cluster & !is_bad_qc)
  removed_by_qc <- sum(is_bad_qc & !is_bad_cluster)
  removed_by_both <- sum(is_bad_cluster & is_bad_qc)

  cat("Filtering with cluster removal:\n")
  cat("- Original cell count:", original_count, "\n")
  cat("- Clusters detected:",
      paste(
            levels(sce[[paste0("clusters_louvain_K", knn)]]),
            collapse = ", "), "\n")
  cat("- Cells remaining:", ncol(sce_filtered), "\n")
  cat(
    "- Clusters removed:",
    paste(clusters_to_remove, collapse = ", "),
    "\n"
  )
  cat(
    "- Removed",
    removed_by_cluster,
    "cells by cluster criterion only\n"
  )
  cat("- Removed", removed_by_qc, "cells by QC criteria only\n")
  cat("- Removed", removed_by_both, "cells by both criteria\n")
  cat(
    "- Total removed:",
    original_count - ncol(sce_filtered),
    "cells (",
    round((original_count - ncol(sce_filtered)) / original_count * 100, 2),
    "%)\n"
  )
} else {
  keep_cells <- !is_bad_qc

  sce_filtered <- sce[, keep_cells]

  cat(
    "No clusters exceed the threshold of",
    cluster_removal_threshold,
    "% bad cells\n"
  )
  cat("Filtering by QC criteria only:\n")
  cat("- Original cell count:", original_count, "\n")
  cat("- Cells remaining:", ncol(sce_filtered), "\n")
  cat(
    "- Removed",
    original_count - ncol(sce_filtered),
    "cells (",
    round((original_count - ncol(sce_filtered)) / original_count * 100, 2),
    "%)\n"
  )
}

# Export outputs
## Unfiltered SCE - store with all the data
saveRDS(sce, file = output_sce_unfiltered)
## Filtered SCE - clean from unrelevant data and store
logcounts(sce_filtered) <- NULL
reducedDim(sce_filtered, type = "PCA") <- NULL
reducedDim(sce_filtered, type = "UMAP") <- NULL
sce_filtered$clusters_louvain_K30 <- NULL
sce_filtered$scDblFinder.class <- NULL
sce_filtered$scDblFinder.score <- NULL
sce_filtered$scDblFinder.weighted <- NULL
sce_filtered$scDblFinder.cxds_score <- NULL
sce_filtered$subsets_Mito_filter_out <- NULL
sce_filtered$discard <- NULL
sce_filtered$low_lib_size <- NULL
sce_filtered$low_n_features <- NULL
sce_filtered$subsets_Mito_detected <- NULL
sce_filtered$subsets_Mito_sum <- NULL
sce_filtered$high_subsets_Mito_percent <- NULL
sce_filtered$sizeFactor <- NULL
rowData(sce_filtered)[["subsets_Mito"]] <- NULL
rowData(sce_filtered)[["scDblFinder.selected"]] <- NULL
metadata(sce_filtered)[["scDblFinder.threshold"]] <- NULL

saveRDS(sce_filtered, file = output_sce_filtered)

## Export plots
ggsave(
  filename = output_qc_plot_transcripts,
  plot = plot_qc_transcripts,
  bg = "white",
  dpi = 300,
  units = "px",
  width = 2228,
  height = 1812
)
ggsave(
  filename = output_qc_plot_doublets,
  plot = plot_doublets,
  bg = "white",
  dpi = 300,
  units = "px",
  width = 2228,
  height = 1812
)
ggsave(
  filename = output_qc_plot_umaps,
  plot = plot_qc_umaps,
  bg = "white",
  dpi = 300,
  units = "px",
  width = 2228,
  height = 1812
)

# Close log connection
sink(type = "message")
sink(type = "output")
close(log)
