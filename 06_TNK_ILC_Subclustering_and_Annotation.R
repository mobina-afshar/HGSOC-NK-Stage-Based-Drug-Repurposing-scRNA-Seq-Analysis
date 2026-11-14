
############################################################
# GSE184880 - HGSOC scRNA-seq Analysis
# Script 06: T/NK/ILC Subclustering and Annotation
#
# Purpose:
# Select T/NK/ILC cells from the CellTypist-annotated immune
# subset, perform subclustering, identify cluster markers,
# and assign the final cell-type annotations used in
# the stage-associated downstream workflow.
#
# Input:
# - config.yaml
# - results/rds/05_HGSOC_immune_CellTypist.rds
#
# Output:
# - Available CellTypist labels
# - T/NK/ILC cluster markers
# - T/NK/ILC cluster sizes
# - Final annotation summary
# - Final T/NK/ILC UMAP
# - results/rds/06_HGSOC_TNK_ILC_Annotated_All.rds
# - results/rds/07_HGSOC_TNK_ILC_Final_Annotated.rds
#
# IMPORTANT:
# - T/NK/ILC selection is based on the original CellTypist labels.
# - Clustering parameters are retrieved from config.yaml.
# - Final cluster annotations reproduce the established
#   annotation worksheet from the original analysis.
# - Clusters 8 and 9 are classified as contaminants and excluded
#   from the final retained object.
# - Cluster 0 and cluster 11 are retained as the main populations
#   used in the downstream stage-associated analysis.
############################################################


############################################################
# 1. Load Project Configuration
############################################################

# Ensure that the script is launched from the repository root.

if (!file.exists("config.yaml")) {
  
  stop(
    "config.yaml was not found. ",
    "Run this script from the project root."
  )
}


# Check and load the YAML package.

if (!requireNamespace("yaml", quietly = TRUE)) {
  
  stop(
    "Package 'yaml' is required. ",
    "Install it using install.packages('yaml')."
  )
}


# Read project configuration.

cfg <- yaml::read_yaml("config.yaml")


# Define project directories.

PROJECT_ROOT <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

RESULTS_DIR <- file.path(
  PROJECT_ROOT,
  cfg$paths$results
)

LOGS_DIR <- file.path(
  PROJECT_ROOT,
  cfg$paths$logs
)

RAW_DIR <- file.path(
  PROJECT_ROOT,
  cfg$paths$raw_data
)

RDS_DIR <- file.path(
  RESULTS_DIR,
  "rds"
)


# Create required output directories.

dir.create(
  RESULTS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  LOGS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  RDS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Set project-wide random seed.

set.seed(cfg$project$seed)


############################################################
# 2. Define Shared Helper Functions
############################################################

# Generate output file paths and create parent directories.

out_dir <- function(...) {
  
  path <- file.path(
    RESULTS_DIR,
    ...
  )
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  return(path)
}


# Check whether required R packages are installed.

require_pkgs <- function(pkgs) {
  
  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  
  if (length(missing) > 0) {
    
    stop(
      "Install missing packages before continuing: ",
      paste(missing, collapse = ", ")
    )
  }
}


# Save R session information for reproducibility.

save_session <- function(tag) {
  
  capture.output(
    sessionInfo(),
    file = file.path(
      LOGS_DIR,
      paste0(tag, "_sessionInfo.txt")
    )
  )
}


############################################################
# 3. Load Required Packages
############################################################

require_pkgs(
  c(
    "Seurat",
    "dplyr",
    "ggplot2"
  )
)

library(Seurat)
library(dplyr)
library(ggplot2)


############################################################
# 4. Load CellTypist-Annotated Immune-Cell Object
############################################################

# Load the immune-cell subset annotated by CellTypist
# in Script 05.

input_rds <- out_dir(
  "rds",
  "05_HGSOC_immune_CellTypist.rds"
)

if (!file.exists(input_rds)) {
  
  stop(
    "Run Script 05 first. Missing input file: ",
    input_rds
  )
}


immune <- readRDS(input_rds)


############################################################
# 5. Identify T/NK/ILC Cells Using CellTypist Labels
############################################################

# Define the exact CellTypist labels used in the original
# source analysis.
#
# IMPORTANT:
# Do not modify these labels without reviewing the actual
# CellTypist output and documenting the reason.

tnk_labels <- c(
  "T cells",
  "T cell",
  "NK cells",
  "NK cell",
  "ILC"
)


# Export all available CellTypist labels for auditing.

available_labels <- sort(
  unique(
    as.character(immune$CellTypist_label)
  )
)


writeLines(
  
  available_labels,
  
  out_dir(
    "tnk_ilc",
    "01_available_celltypist_labels.txt"
  )
)


# Select cells matching the predefined T/NK/ILC labels.

selected_cells <- colnames(immune)[
  
  as.character(immune$CellTypist_label) %in% tnk_labels
]


# Stop if no cells match the expected labels.

if (length(selected_cells) == 0) {
  
  stop(
    "No cells matched the source CellTypist T/NK/ILC labels. ",
    "Review results/tnk_ilc/01_available_celltypist_labels.txt; ",
    "do not guess label mappings."
  )
}


# Create the T/NK/ILC subset.

tnk <- subset(
  
  immune,
  
  cells = selected_cells
)


############################################################
# 6. Prepare T/NK/ILC Expression Data
############################################################

# Set the RNA assay as the default assay.

DefaultAssay(tnk) <- "RNA"


# Normalize expression data using LogNormalize.

tnk <- NormalizeData(
  
  tnk,
  
  normalization.method = "LogNormalize",
  
  scale.factor = 10000,
  
  verbose = FALSE
)


# Identify highly variable genes.

tnk <- FindVariableFeatures(
  
  tnk,
  
  selection.method = "vst",
  
  nfeatures = cfg$seurat$n_hvg,
  
  verbose = FALSE
)


# Scale highly variable genes.

tnk <- ScaleData(
  
  tnk,
  
  features = VariableFeatures(tnk),
  
  verbose = FALSE
)


############################################################
# 7. Perform PCA
############################################################

# Run PCA using the highly variable genes.
#
# The number of computed principal components is set to 30.
# The number used in downstream analysis is configured
# separately through cfg$seurat$n_pcs_tnk_ilc.

tnk <- RunPCA(
  
  tnk,
  
  features = VariableFeatures(tnk),
  
  npcs = 30,
  
  seed.use = cfg$project$seed,
  
  verbose = FALSE
)


############################################################
# 8. Perform T/NK/ILC Subclustering
############################################################

# Construct the nearest-neighbor graph using the configured
# number of principal components.

tnk <- FindNeighbors(
  
  tnk,
  
  reduction = "pca",
  
  dims = seq_len(cfg$seurat$n_pcs_tnk_ilc),
  
  verbose = FALSE
)


# Identify cell clusters using the configured resolution.

tnk <- FindClusters(
  
  tnk,
  
  resolution = cfg$seurat$tnk_ilc_resolution,
  
  algorithm = 1,
  
  random.seed = cfg$project$seed,
  
  verbose = FALSE
)


# Generate UMAP coordinates for visualization.

tnk <- RunUMAP(
  
  tnk,
  
  reduction = "pca",
  
  dims = seq_len(cfg$seurat$n_pcs_tnk_ilc),
  
  seed.use = cfg$project$seed,
  
  verbose = FALSE
)


############################################################
# 9. Identify Cluster Marker Genes
############################################################

# Identify positively enriched marker genes for each
# T/NK/ILC subcluster.
#
# These markers provide an auditable basis for reviewing
# the biological interpretation of the clusters.

markers <- FindAllMarkers(
  
  tnk,
  
  assay = "RNA",
  
  only.pos = TRUE,
  
  min.pct = 0.25,
  
  logfc.threshold = 0.25,
  
  test.use = "wilcox",
  
  random.seed = cfg$project$seed
)


# Export cluster marker results.

write.csv(
  
  markers,
  
  out_dir(
    "tnk_ilc",
    "02_tnk_ilc_cluster_markers.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 10. Summarize T/NK/ILC Cluster Sizes
############################################################

# Calculate the number of cells assigned to each subcluster.

cluster_sizes <- as.data.frame(
  table(tnk$seurat_clusters)
)


# Export cluster sizes.

write.csv(
  
  cluster_sizes,
  
  out_dir(
    "tnk_ilc",
    "03_tnk_ilc_cluster_sizes.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 11. Assign Final Cluster Annotations
############################################################

# Final cluster labels reproduce the established annotation
# worksheet from the original analysis.
#
# IMPORTANT:
# Cluster identities are linked to the established clustering
# result. Do not change the mapping without reviewing marker
# expression and confirming cluster correspondence.

final_annotation <- c(
  
  "0" = "NK/ILC3-like",
  
  "1" = "Activated cytotoxic T/NK-like",
  
  "2" = "Activated/Exhausted NK/T-like",
  
  "3" = "NK/ILC-like",
  
  "4" = "Cycling NK/ILC",
  
  "5" = "Mature cytotoxic NK",
  
  "6" = "IFN-stimulated NK/ILC",
  
  "7" = "Cycling NK/ILC",
  
  "8" = "Stromal/Epithelial contaminant",
  
  "9" = "Myeloid contaminant",
  
  "10" = "ILC-like / ILC1-like",
  
  "11" = "ILC3"
)


# Identify clusters present in the current dataset.

clusters_present <- unique(
  as.character(tnk$seurat_clusters)
)


# Check whether every detected cluster has a predefined label.

if (!all(clusters_present %in% names(final_annotation))) {
  
  stop(
    "Unexpected T/NK/ILC clusters detected. ",
    "Review cluster markers and update annotation explicitly."
  )
}


# Assign final annotations to individual cells.

tnk$Final_Annotation <- unname(
  
  final_annotation[
    as.character(tnk$seurat_clusters)
  ]
)


############################################################
# 12. Assign Annotation Confidence
############################################################

# Define confidence levels for the established cluster labels.
#
# These confidence levels reproduce the original annotation
# worksheet and are not independently calculated by Seurat.

confidence_map <- c(
  
  "0" = "Moderate",
  
  "1" = "Moderate",
  
  "2" = "Low-Moderate",
  
  "3" = "Low-Moderate",
  
  "4" = "High",
  
  "5" = "Very High",
  
  "6" = "High for functional state",
  
  "7" = "Very High",
  
  "8" = "Very High",
  
  "9" = "Very High",
  
  "10" = "Moderate",
  
  "11" = "Very High"
)


# Assign confidence values to individual cells
# according to their cluster identity.

confidence_values <- unname(
  
  confidence_map[
    as.character(tnk$seurat_clusters)
  ]
)


# Check for missing confidence assignments.

if (anyNA(confidence_values)) {
  
  stop(
    "Missing confidence annotations for some clusters."
  )
}


# Add annotation confidence to Seurat metadata.

tnk$Annotation_Confidence <- confidence_values


############################################################
# 13. Remove Contaminant Clusters
############################################################

# Define clusters classified as stromal/epithelial
# and myeloid contaminants.

contaminant_clusters <- c(
  "8",
  "9"
)


# Retain all other annotated T/NK/ILC clusters.

tnk_retained <- subset(
  
  tnk,
  
  subset = !seurat_clusters %in% contaminant_clusters
)


# Set final biological annotations as active identities.

Idents(tnk_retained) <- "Final_Annotation"


############################################################
# 14. Generate Final Annotation Summary
############################################################

# Summarize the number of cells in each cluster,
# final annotation, and annotation-confidence category.
#
# The summary is generated from the complete T/NK/ILC object,
# before contaminant removal.

annotation_summary <- tnk@meta.data |>
  
  dplyr::count(
    
    seurat_clusters,
    
    Final_Annotation,
    
    Annotation_Confidence,
    
    name = "cells"
  )


# Export the final annotation summary.

write.csv(
  
  annotation_summary,
  
  out_dir(
    "tnk_ilc",
    "04_final_tnk_ilc_annotation_summary.csv"
  ),
  
  row.names = FALSE
)


# Report the total number of annotated cells before
# contaminant removal.

cat(
  "Total annotated T/NK/ILC cells:",
  sum(annotation_summary$cells),
  "\n"
)

cat(
  "Total cells in T/NK/ILC object:",
  ncol(tnk),
  "\n"
)


############################################################
# 15. Visualize Final T/NK/ILC Annotations
############################################################

# Generate UMAP visualization using the retained
# T/NK/ILC populations.
#
# Clusters 8 and 9 are excluded from this visualization.

p <- DimPlot(
  
  tnk_retained,
  
  reduction = "umap",
  
  group.by = "Final_Annotation",
  
  label = TRUE,
  
  repel = TRUE
) +
  
  theme_classic()


# Save final T/NK/ILC UMAP.

ggsave(
  
  out_dir(
    "tnk_ilc",
    "05_final_tnk_ilc_umap.png"
  ),
  
  p,
  
  width = 14,
  
  height = 9,
  
  dpi = 300
)


############################################################
# 16. Save Final Seurat Objects
############################################################

# Save the complete annotated T/NK/ILC object,
# including contaminant clusters.

saveRDS(
  
  tnk,
  
  out_dir(
    "rds",
    "06_HGSOC_TNK_ILC_Annotated_All.rds"
  )
)


# Save the retained T/NK/ILC object after removing
# contaminant clusters 8 and 9.

saveRDS(
  
  tnk_retained,
  
  out_dir(
    "rds",
    "07_HGSOC_TNK_ILC_Final_Annotated.rds"
  )
)


############################################################
# 17. Report Final Results
############################################################

cat("\n")

cat("============================================================\n")
cat("SCRIPT 06: T/NK/ILC SUBCLUSTERING COMPLETED\n")
cat("============================================================\n")

cat(
  "Initial immune cells:",
  ncol(immune),
  "\n"
)

cat(
  "Selected T/NK/ILC cells:",
  ncol(tnk),
  "\n"
)

cat(
  "Number of subclusters:",
  length(unique(tnk$seurat_clusters)),
  "\n"
)

cat(
  "Retained cells after contaminant removal:",
  ncol(tnk_retained),
  "\n"
)

cat(
  "Removed contaminant clusters:",
  paste(contaminant_clusters, collapse = ", "),
  "\n"
)

cat(
  "Complete annotated object: ",
  out_dir(
    "rds",
    "06_HGSOC_TNK_ILC_Annotated_All.rds"
  ),
  "\n",
  sep = ""
)

cat(
  "Final retained object: ",
  out_dir(
    "rds",
    "07_HGSOC_TNK_ILC_Final_Annotated.rds"
  ),
  "\n",
  sep = ""
)


############################################################
# 18. Save Session Information
############################################################

save_session("06_tnk_ilc")


############################################################
# END OF SCRIPT 06
############################################################
