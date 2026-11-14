
############################################################
# GSE184880 - HGSOC scRNA-seq Analysis
# Script 08: Differential Expression Analysis with DESeq2
#
# Workflow: Patient-Level Differential Expression Analysis
#
# Purpose:
# Perform population-specific differential expression analysis
# using patient-level pseudobulk raw counts.
#
# Populations:
# - Cluster 0  -> NK-like
# - Cluster 11 -> ILC3
#
# Comparison:
# - Late vs Early
# - Reference group: Early
#
# Input:
# - config.yaml
# - results/rds/09_pseudobulk_counts.rds
# - results/rds/10_pseudobulk_metadata.rds
#
# Output:
# - Population-specific DESeq2 results
# - Population-specific normalized counts
# - Population-specific DESeq2 objects
# - Combined DESeq2 results
# - DESeq2 analysis summary
# - DESeq2 analysis parameters
# - results/rds/DESeq2_Results_All_Populations.rds
# - logs/08_DESeq2_sessionInfo.txt
#
# IMPORTANT:
# - The unit of analysis is the patient-level pseudobulk sample.
# - NK-like and ILC3 are analyzed separately.
# - Raw integer pseudobulk counts are used as input.
# - Early is the reference group.
# - The gene filter is rowSums(counts >= 10) >= 2.
# - No log2FC-based gene filtering is applied in this script.
# - Positive log2FC indicates higher expression in Late.
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

RDS_DIR <- file.path(
  RESULTS_DIR,
  "rds"
)

DE_DIR <- file.path(
  RESULTS_DIR,
  "differential_expression"
)


# Create required output directories.

dir.create(
  DE_DIR,
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
# 2. Define Analysis Parameters
############################################################

# Define the populations included in Workflow 2.

populations <- c(
  "NK-like",
  "ILC3"
)


# Define the reference and comparison groups.

reference_group <- "Early"

comparison_group <- "Late"


# Define the adjusted p-value threshold.

padj_cutoff <- 0.05


# Define gene-level count filtering parameters.
#
# Retain genes with at least 10 counts in at least
# 2 pseudobulk samples within each population.

min_count <- 10

min_samples <- 2


############################################################
# 3. Define Shared Helper Functions
############################################################

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


# Generate output file paths and create parent directories.

output_file <- function(...) {
  
  path <- file.path(
    DE_DIR,
    ...
  )
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  return(path)
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
# 4. Load Required Packages
############################################################

require_pkgs(
  c(
    "DESeq2",
    "Matrix"
  )
)

library(DESeq2)
library(Matrix)


############################################################
# 5. Load Pseudobulk Data
############################################################

# Define input file paths.

counts_file <- file.path(
  RDS_DIR,
  "09_pseudobulk_counts.rds"
)

metadata_file <- file.path(
  RDS_DIR,
  "10_pseudobulk_metadata.rds"
)


# Check whether the pseudobulk count matrix exists.

if (!file.exists(counts_file)) {
  
  stop(
    "Pseudobulk count matrix not found. ",
    "Run Script 07 first."
  )
}


# Check whether pseudobulk metadata exists.

if (!file.exists(metadata_file)) {
  
  stop(
    "Pseudobulk metadata not found. ",
    "Run Script 07 first."
  )
}


# Load pseudobulk count matrix and metadata.

pb_counts <- readRDS(counts_file)

pb_meta <- readRDS(metadata_file)


############################################################
# 6. Validate Pseudobulk Input Data
############################################################

# Check count matrix structure.

if (
  !inherits(pb_counts, "Matrix") &&
  !is.matrix(pb_counts)
) {
  
  stop(
    "Pseudobulk counts must be a matrix."
  )
}


# Check metadata structure.

if (!is.data.frame(pb_meta)) {
  
  stop(
    "Pseudobulk metadata must be a data frame."
  )
}


# Define required metadata columns.

required_metadata <- c(
  "pseudobulk_id",
  "Population",
  "sample",
  "stage_group"
)


# Identify missing metadata columns.

missing_metadata <- setdiff(
  required_metadata,
  colnames(pb_meta)
)


if (length(missing_metadata) > 0) {
  
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}


# Verify that count matrix identifiers are available.

if (
  is.null(rownames(pb_counts)) ||
  is.null(colnames(pb_counts))
) {
  
  stop(
    "Count matrix must have gene and sample identifiers."
  )
}


# Check for duplicated gene identifiers.

if (anyDuplicated(rownames(pb_counts))) {
  
  stop(
    "Duplicated gene identifiers detected."
  )
}


# Check for duplicated pseudobulk identifiers.

if (anyDuplicated(colnames(pb_counts))) {
  
  stop(
    "Duplicated pseudobulk identifiers detected."
  )
}


# Verify consistency between count matrix and metadata identifiers.

if (!setequal(
  colnames(pb_counts),
  rownames(pb_meta)
)) {
  
  stop(
    "Count matrix columns and metadata row names do not match."
  )
}


# Reorder metadata to match the count matrix.

pb_meta <- pb_meta[
  colnames(pb_counts),
  ,
  drop = FALSE
]


# Confirm identifier consistency.

stopifnot(
  identical(
    rownames(pb_meta),
    colnames(pb_counts)
  )
)


# Check for missing count values.

if (anyNA(pb_counts)) {
  
  stop(
    "Missing values detected in pseudobulk counts."
  )
}


# Check for negative count values.

if (any(pb_counts < 0)) {
  
  stop(
    "Negative values detected in pseudobulk counts."
  )
}


# Verify that all count values are integers.

if (any(
  abs(pb_counts - round(pb_counts)) > 1e-8
)) {
  
  stop(
    "Pseudobulk counts must be integer-valued raw counts."
  )
}


# Check for missing clinical stage information.

if (anyNA(pb_meta$stage_group)) {
  
  stop(
    "Missing stage_group values detected."
  )
}


# Verify the expected clinical stage groups.

if (!all(
  as.character(pb_meta$stage_group) %in%
  c("Early", "Late")
)) {
  
  stop(
    "Unexpected stage_group values detected."
  )
}


# Display input data summary.

cat("\n")
cat("============================================================\n")
cat("INPUT DATA SUMMARY\n")
cat("============================================================\n")

cat(
  "Number of genes:",
  nrow(pb_counts),
  "\n"
)

cat(
  "Number of pseudobulk samples:",
  ncol(pb_counts),
  "\n"
)

cat("\nPopulation distribution:\n")

print(
  table(pb_meta$Population)
)

cat("\nStage distribution:\n")

print(
  table(pb_meta$stage_group)
)


############################################################
# 7. Define Differential Expression Analysis Function
############################################################

run_deseq2 <- function(population_name) {
  
  cat("\n")
  
  cat("============================================================\n")
  
  cat(
    "Population:",
    population_name,
    "\n"
  )
  
  cat("============================================================\n")
  
  
  ##########################################################
  # 7.1 Select Population-Specific Pseudobulk Samples
  ##########################################################
  
  # Select metadata records corresponding to the current population.
  
  selected_meta <- pb_meta[
    pb_meta$Population == population_name,
    ,
    drop = FALSE
  ]
  
  
  # Skip the population if no pseudobulk samples are available.
  
  if (nrow(selected_meta) == 0) {
    
    message(
      "No pseudobulk samples available for ",
      population_name,
      ". Skipping."
    )
    
    return(NULL)
  }
  
  
  # Extract the corresponding pseudobulk count matrix.
  
  selected_counts <- pb_counts[
    ,
    rownames(selected_meta),
    drop = FALSE
  ]
  
  
  ##########################################################
  # 7.2 Check Biological Replicates
  ##########################################################
  
  # Set Early as the reference group and Late as the comparison group.
  
  selected_meta$stage_group <- factor(
    
    as.character(selected_meta$stage_group),
    
    levels = c(
      reference_group,
      comparison_group
    )
  )
  
  
  # Count biological replicates in each stage group.
  
  replicate_table <- table(
    selected_meta$stage_group
  )
  
  
  cat("\nBiological replicates:\n")
  
  print(replicate_table)
  
  
  # Both stage groups must be present for differential
  # expression analysis.
  
  if (any(replicate_table == 0)) {
    
    message(
      "Both Early and Late groups are required for ",
      population_name,
      ". Differential expression skipped."
    )
    
    return(NULL)
  }
  
  
  # Warn when either group contains fewer than two
  # biological replicates.
  
  if (any(replicate_table < 2)) {
    
    warning(
      "Fewer than two biological replicates in at least one group for ",
      population_name,
      ". Results will have limited statistical reliability."
    )
  }
  
  
  ##########################################################
  # 7.3 Filter Low-Count Genes
  ##########################################################
  
  # Retain genes with at least min_count counts
  # in at least min_samples pseudobulk samples.
  
  keep_genes <- Matrix::rowSums(
    selected_counts >= min_count
  ) >= min_samples
  
  
  # Apply gene-level filtering.
  
  count_sub <- selected_counts[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  
  # Report gene numbers before and after filtering.
  
  cat(
    "\nGenes before filtering:",
    nrow(selected_counts),
    "\n"
  )
  
  cat(
    "Genes after filtering:",
    nrow(count_sub),
    "\n"
  )
  
  
  # Check whether enough genes remain for analysis.
  
  if (nrow(count_sub) < 2) {
    
    warning(
      "Insufficient genes after filtering for ",
      population_name,
      ". Differential expression skipped."
    )
    
    return(NULL)
  }
  
  
  ##########################################################
  # 7.4 Prepare DESeq2 Input
  ##########################################################
  
  # Convert the filtered pseudobulk count matrix to
  # a standard integer matrix.
  
  count_matrix <- as.matrix(
    count_sub
  )
  
  storage.mode(count_matrix) <- "integer"
  
  
  # Prepare sample-level metadata.
  
  col_data <- as.data.frame(
    selected_meta
  )
  
  
  # Set metadata row names to pseudobulk identifiers.
  
  rownames(col_data) <- rownames(
    selected_meta
  )
  
  
  # Reorder metadata to match count matrix columns.
  
  col_data <- col_data[
    colnames(count_matrix),
    ,
    drop = FALSE
  ]
  
  
  # Verify sample order consistency.
  
  stopifnot(
    identical(
      rownames(col_data),
      colnames(count_matrix)
    )
  )
  
  
  # Drop unused stage factor levels.
  
  col_data$stage_group <- droplevels(
    col_data$stage_group
  )
  
  
  ##########################################################
  # 7.5 Construct DESeq2 Dataset
  ##########################################################
  
  # Create the DESeq2 dataset using patient-level pseudobulk counts.
  #
  # Design:
  # ~ stage_group
  #
  # Early is the reference level.
  
  dds <- DESeqDataSetFromMatrix(
    
    countData = count_matrix,
    
    colData = col_data,
    
    design = ~ stage_group
  )
  
  
  ##########################################################
  # 7.6 Run Differential Expression Analysis
  ##########################################################
  
  # Estimate size factors, dispersions and fit the
  # negative binomial model.
  
  dds <- DESeq(dds)
  
  
  ##########################################################
  # 7.7 Extract Late vs Early Results
  ##########################################################
  
  # Extract the differential expression contrast.
  #
  # Positive log2FoldChange:
  # Higher expression in Late.
  #
  # Negative log2FoldChange:
  # Higher expression in Early.
  
  res <- results(
    
    dds,
    
    contrast = c(
      "stage_group",
      comparison_group,
      reference_group
    ),
    
    alpha = padj_cutoff
  )
  
  
  ##########################################################
  # 7.8 Format Differential Expression Results
  ##########################################################
  
  # Convert DESeq2 results to a data frame.
  
  res_df <- as.data.frame(
    res
  )
  
  
  # Add gene identifiers.
  
  res_df$gene <- rownames(
    res_df
  )
  
  
  # Retain the required result columns.
  
  res_df <- res_df[
    ,
    c(
      "gene",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )
  ]
  
  
  # Order results by adjusted p-value.
  # Missing adjusted p-values are placed at the end.
  
  res_df <- res_df[
    order(
      is.na(res_df$padj),
      res_df$padj
    ),
    ,
    drop = FALSE
  ]
  
  
  ##########################################################
  # 7.9 Add Population Information
  ##########################################################
  
  # Add the population identifier to each result.
  
  res_df$Population <- population_name
  
  
  ##########################################################
  # 7.10 Save Population-Specific DESeq2 Results
  ##########################################################
  
  # Generate a file-safe population name.
  
  population_tag <- gsub(
    "-",
    "_",
    population_name
  )
  
  
  # Export the complete DESeq2 result table.
  
  write.csv(
    
    res_df,
    
    output_file(
      paste0(
        population_tag,
        "_DESeq2_all_results.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  ##########################################################
  # 7.11 Save Normalized Counts
  ##########################################################
  
  # Extract DESeq2 size-factor normalized counts.
  
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  
  # Export normalized counts.
  
  write.csv(
    
    as.data.frame(normalized_counts),
    
    output_file(
      paste0(
        population_tag,
        "_normalized_counts.csv"
      )
    )
  )
  
  
  ##########################################################
  # 7.12 Save Population-Specific DESeq2 Objects
  ##########################################################
  
  # Save the DESeq2 dataset object.
  
  saveRDS(
    
    dds,
    
    output_file(
      paste0(
        population_tag,
        "_DESeq2_object.rds"
      )
    )
  )
  
  
  # Save the DESeq2 results object.
  
  saveRDS(
    
    res,
    
    output_file(
      paste0(
        population_tag,
        "_DESeq2_results.rds"
      )
    )
  )
  
  
  ##########################################################
  # 7.13 Generate Population-Specific Analysis Summary
  ##########################################################
  
  # Summarize the number of biological replicates
  # and genes retained after filtering.
  
  summary_df <- data.frame(
    
    Population = population_name,
    
    Early_replicates = as.integer(
      replicate_table[reference_group]
    ),
    
    Late_replicates = as.integer(
      replicate_table[comparison_group]
    ),
    
    Genes_before_filtering = nrow(selected_counts),
    
    Genes_after_filtering = nrow(count_sub),
    
    stringsAsFactors = FALSE
  )
  
  
  cat("\n")
  
  cat("===== DESeq2 SUMMARY =====\n")
  
  print(summary_df)
  
  
  ##########################################################
  # 7.14 Return Population-Specific Results
  ##########################################################
  
  return(
    list(
      
      results = res_df,
      
      summary = summary_df,
      
      dds = dds
    )
  )
}


############################################################
# 8. Run DESeq2 for Each Population
############################################################

# Perform separate differential expression analyses
# for NK-like and ILC3 populations.

deseq_results <- lapply(
  
  populations,
  
  run_deseq2
)


# Preserve population names in the result list.

names(deseq_results) <- populations


############################################################
# 9. Combine Successful Population Results
############################################################

# Identify populations for which DESeq2 analysis completed.

successful_results <- Filter(
  
  Negate(is.null),
  
  deseq_results
)


# Stop if no population could be analyzed.

if (length(successful_results) == 0) {
  
  stop(
    "No population could be analyzed. ",
    "Check the Early/Late sample distribution."
  )
}


# Combine population-specific differential expression tables.

all_results <- do.call(
  
  rbind,
  
  lapply(
    successful_results,
    function(x) x$results
  )
)


# Combine population-specific analysis summaries.

all_summaries <- do.call(
  
  rbind,
  
  lapply(
    successful_results,
    function(x) x$summary
  )
)


############################################################
# 10. Save Combined Differential Expression Results
############################################################

# Export combined differential expression results.

write.csv(
  
  all_results,
  
  output_file(
    "All_Populations_DESeq2_results.csv"
  ),
  
  row.names = FALSE
)


# Export combined analysis summary.

write.csv(
  
  all_summaries,
  
  output_file(
    "DESeq2_Analysis_Summary.csv"
  ),
  
  row.names = FALSE
)


# Save the complete named result list.
#
# Each population contains:
# - $results: differential expression result table
# - $summary: population-specific analysis summary
# - $dds: DESeq2 dataset object
#
# This structure is required by downstream Scripts 09 and 10.

saveRDS(
  
  deseq_results,
  
  output_file(
    "DESeq2_Results_All_Populations.rds"
  )
)


############################################################
# 11. Save Analysis Parameters
############################################################

# Record the main parameters used in the analysis.

analysis_parameters <- list(
  
  dataset = "GSE184880",
  
  method = "DESeq2",
  
  unit_of_analysis = "Patient-level pseudobulk",
  
  design = "~ stage_group",
  
  contrast = c(
    "stage_group",
    "Late",
    "Early"
  ),
  
  reference_group = reference_group,
  
  comparison_group = comparison_group,
  
  gene_filter = "rowSums(counts >= 10) >= 2",
  
  min_count = min_count,
  
  min_samples = min_samples,
  
  padj_cutoff = padj_cutoff,
  
  populations = populations
)


# Save analysis parameters.

saveRDS(
  
  analysis_parameters,
  
  output_file(
    "DESeq2_Analysis_Parameters.rds"
  )
)


############################################################
# 12. Save Session Information
############################################################

save_session("08_DESeq2")


############################################################
# 13. Final Analysis Summary
############################################################

cat("\n")

cat("============================================================\n")

cat("SCRIPT 08: DESeq2 ANALYSIS COMPLETED\n")

cat("============================================================\n")


# Display the population-specific analysis summary.

cat("\nDESeq2 summary:\n")

print(all_summaries)


# Display the output directory.

cat("\nOutput directory:\n")

cat(
  DE_DIR,
  "\n"
)


# Identify the next analysis step.

cat(
  "\nNext step: Script 09 - Hallmark GSEA\n"
)


############################################################
# END OF SCRIPT 08
############################################################
