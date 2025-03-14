threads <- snakemake@params[["threads"]]
sample_name <- snakemake@params[["sample"]]
counts_path <- snakemake@params[["counts_dir"]]
vdj_b_path <- snakemake@params[["vdj_b"]]
vdj_t_path <- snakemake@params[["vdj_t"]]
output_sce <- snakemake@output[["sce"]]

# Parallelization
library(BiocParallel)
bpp <- MulticoreParam(workers = threads)

# Logging
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load the counts
library(DropletUtils)
sce <- read10xCounts(
  counts_path,
  sample.names = sample_name,
  col.names = TRUE,
  BPPARAM = bpp
)
# Barcode unification
sce$Sample_Barcode <- paste(sce$Sample, sce$Barcode, sep = "_")
colnames(sce) <- sce$Sample_Barcode

# Load V(D)J BCR data
library(scRepertoire)
contigs_bcr <- read.csv(vdj_b_path)
clonotypes_bcr <- combineBCR(
  contigs_bcr,
  samples = sample_name,
  call.related.clones = TRUE,
  threshold = 0.85,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = TRUE,
  filterNonproductive = TRUE
)

# Load V(D)J TCR data
contigs_tcr <- read.csv(vdj_t_path)
clonotypes_tcr <- combineTCR(
  contigs_tcr,
  samples = sample_name,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = TRUE,
  filterNonproductive = TRUE
)

# Combine counts and V(D)J data
sce <- combineExpression(
  c(clonotypes_bcr, clonotypes_tcr),
  sce,
  cloneCall = "strict",
  chain = "both",
  group.by = NULL,
  proportion = TRUE,
  filterNA = FALSE,
  cloneSize = c(
    Rare = 1e-04,
    Small = 0.001,
    Medium = 0.01,
    Large = 0.1,
    Hyperexpanded = 1
  )
)

# Log output
print("Resulting SCE object:")
print(sce)
print("colData names in the SCE object:")
print(colnames(colData(sce)))
print(paste("Saving the SCE as", output_sce))
# Save SingleCellExperiment object as RDS
saveRDS(sce, file = output_sce)

# Close log connection
sink(type = "message")
sink(type = "output")
close(log)