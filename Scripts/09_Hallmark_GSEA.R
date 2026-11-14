
############################################################
# GSE184880 - HGSOC scRNA-seq Analysis
# Script 09: Preranked Hallmark GSEA
#
# Workflow: Stage-Based Drug Repurposing in HGSOC
#
# Purpose:
# Perform preranked Hallmark Gene Set Enrichment Analysis
# using patient-level pseudobulk DESeq2 results.
#
# Cell populations:
# - Cluster 0  -> NK-like
# - Cluster 11 -> ILC3
#
# Comparison:
# - Late vs Early
#
# Main downstream population:
# - NK-like -> CMap / Drug Repurposing
#
# ILC3:
# - Complementary pathway analysis only
#
# Input:
# - config.yaml
# - results/rds/DESeq2_Results_All_Populations.rds
#
# Output:
# - Population-specific gene rankings
# - Population-specific Hallmark GSEA results
# - FDR-filtered and exploratory pathway tables
# - Combined Hallmark GSEA results
# - GSEA results and parameters in RDS format
# - logs/09_hallmark_gsea_sessionInfo.txt
#
# IMPORTANT:
# - All genes with valid DESeq2 Wald statistics are ranked.
# - No adjusted P-value filter is applied before GSEA.
# - Positive Wald statistics indicate higher expression in Late.
# - Positive NES indicates enrichment in Late.
# - Negative NES indicates enrichment in Early.
# - NK-like is the primary population for downstream
#   drug-repurposing analysis.
############################################################


############################################################
# 1. Load Project Configuration
############################################################

if (!file.exists("config.yaml")) {
  
  stop(
    "config.yaml was not found. ",
    "Run this script from the project root."
  )
}


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

GSEA_DIR <- file.path(
  RESULTS_DIR,
  "gsea"
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

dir.create(
  GSEA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# Set project-wide random seed.

set.seed(cfg$project$seed)


############################################################
# 2. Define Shared Helper Functions
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
    "fgsea",
    "msigdbr"
  )
)

library(fgsea)
library(msigdbr)


############################################################
# 4. Validate GSEA Parameters
############################################################

# Define the required GSEA parameters in config.yaml.

required_gsea_parameters <- c(
  "min_size",
  "max_size",
  "eps",
  "fdr_cutoff",
  "exploratory_p_cutoff"
)


# Check for missing configuration parameters.

missing_gsea_parameters <- setdiff(
  required_gsea_parameters,
  names(cfg$gsea)
)

if (length(missing_gsea_parameters) > 0) {
  
  stop(
    "Missing GSEA parameters in config.yaml: ",
    paste(
      missing_gsea_parameters,
      collapse = ", "
    )
  )
}


# Load GSEA parameters.

min_size <- as.integer(
  cfg$gsea$min_size
)

max_size <- as.integer(
  cfg$gsea$max_size
)

eps_value <- as.numeric(
  cfg$gsea$eps
)

fdr_cutoff <- as.numeric(
  cfg$gsea$fdr_cutoff
)

exploratory_p_cutoff <- as.numeric(
  cfg$gsea$exploratory_p_cutoff
)


# Validate gene set size parameters.

if (min_size < 1 || max_size <= min_size) {
  
  stop(
    "Invalid GSEA min_size or max_size."
  )
}


# Validate fgsea eps parameter.

if (eps_value < 0) {
  
  stop(
    "GSEA eps must be non-negative."
  )
}


# Validate FDR threshold.

if (fdr_cutoff <= 0 || fdr_cutoff > 1) {
  
  stop(
    "Invalid GSEA FDR cutoff."
  )
}


# Validate exploratory nominal P-value threshold.

if (
  exploratory_p_cutoff <= 0 ||
  exploratory_p_cutoff > 1
) {
  
  stop(
    "Invalid exploratory P-value cutoff."
  )
}


############################################################
# 5. Load DESeq2 Results
############################################################

cat("\n")

cat("============================================================\n")

cat("LOADING DESeq2 RESULTS\n")

cat("============================================================\n")


# Define input file.

de_file <- file.path(
  RDS_DIR,
  "DESeq2_Results_All_Populations.rds"
)


# Check whether DESeq2 results are available.

if (!file.exists(de_file)) {
  
  stop(
    "DESeq2 results not found. ",
    "Run Script 08 first."
  )
}


# Load the complete DESeq2 result list.

de_list <- readRDS(de_file)


# Validate expected population names.

expected_populations <- c(
  "NK-like",
  "ILC3"
)

if (!all(
  expected_populations %in% names(de_list)
)) {
  
  stop(
    "Expected populations NK-like and ILC3 ",
    "were not found in the DESeq2 results."
  )
}


# Extract population-specific DESeq2 result tables.

nk_de <- de_list[["NK-like"]]$results

ilc3_de <- de_list[["ILC3"]]$results


# Check that both result tables are available.

if (
  is.null(nk_de) ||
  is.null(ilc3_de)
) {
  
  stop(
    "DESeq2 result tables are missing ",
    "for one or more populations."
  )
}


# Display the number of genes in each population.

cat(
  "\nNK-like genes:",
  nrow(nk_de),
  "\n"
)

cat(
  "ILC3 genes:",
  nrow(ilc3_de),
  "\n"
)


############################################################
# 6. Validate Gene Identifiers
############################################################

# Define required columns in DESeq2 result tables.

required_columns <- c(
  "gene",
  "stat"
)


# Validate each population separately.

for (population_name in expected_populations) {
  
  tab <- de_list[[population_name]]$results
  
  
  # Check required columns.
  
  missing_columns <- setdiff(
    required_columns,
    colnames(tab)
  )
  
  if (length(missing_columns) > 0) {
    
    stop(
      "Missing required columns in ",
      population_name,
      " DESeq2 results: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  
  # Check for duplicated gene identifiers.
  
  if (anyDuplicated(tab$gene)) {
    
    stop(
      "Duplicated gene identifiers detected in ",
      population_name,
      ". Resolve gene identifiers before GSEA."
    )
  }
  
  
  # Check whether Ensembl gene IDs are being used.
  
  if (all(grepl("^ENSG", tab$gene))) {
    
    stop(
      "Ensembl gene IDs detected in ",
      population_name,
      ". Hallmark GSEA requires gene symbols ",
      "for this workflow. Perform identifier mapping first."
    )
  }
}


############################################################
# 7. Load MSigDB Hallmark Gene Sets
############################################################

cat("\n")

cat("============================================================\n")

cat("LOADING HALLMARK GENE SETS\n")

cat("============================================================\n")


# Retrieve the Hallmark collection for Homo sapiens.

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)


# Check required MSigDB columns.

required_hallmark_columns <- c(
  "gs_name",
  "gene_symbol"
)

if (!all(
  required_hallmark_columns %in%
  colnames(hallmark)
)) {
  
  stop(
    "Required Hallmark columns were not found ",
    "in the msigdbr output."
  )
}


# Remove missing and empty gene symbols.

hallmark <- hallmark[
  !is.na(hallmark$gene_symbol) &
    hallmark$gene_symbol != "" &
    !is.na(hallmark$gs_name),
  ,
  drop = FALSE
]


# Construct the pathway list.

pathways <- split(
  hallmark$gene_symbol,
  hallmark$gs_name
)


# Remove duplicated genes within each pathway.

pathways <- lapply(
  pathways,
  unique
)


# Report the number of Hallmark pathways.

cat(
  "Hallmark pathways loaded:",
  length(pathways),
  "\n"
)


# Stop if no Hallmark pathways were retrieved.

if (length(pathways) == 0) {
  
  stop(
    "No Hallmark pathways were retrieved."
  )
}


############################################################
# 8. Prepare Complete Gene Ranking
############################################################

# Rank all genes using the DESeq2 Wald statistic.
#
# Positive Wald statistic:
# Higher expression in Late compared with Early.
#
# Negative Wald statistic:
# Higher expression in Early compared with Late.
#
# Genes are NOT filtered by adjusted P-value.
# This preserves the complete preranked signature.


make_ranks <- function(tab, population_name) {
  
  cat("\n")
  
  cat(
    "Preparing gene ranking for:",
    population_name,
    "\n"
  )
  
  
  # Extract Wald statistics and assign gene identifiers.
  
  x <- tab$stat
  
  names(x) <- as.character(
    tab$gene
  )
  
  
  # Remove missing and non-finite statistics.
  
  valid <- !is.na(x) &
    is.finite(x) &
    !is.na(names(x)) &
    names(x) != ""
  
  x <- x[valid]
  
  
  # Check for duplicated gene identifiers.
  
  if (anyDuplicated(names(x))) {
    
    stop(
      "Duplicated gene identifiers in ",
      population_name,
      " ranking."
    )
  }
  
  
  # Rank all valid genes in descending order.
  
  ranks <- sort(
    x,
    decreasing = TRUE
  )
  
  
  # Identify genes shared between the ranking
  # and Hallmark gene sets.
  
  pathway_genes <- unique(
    unlist(
      pathways,
      use.names = FALSE
    )
  )
  
  overlapping_genes <- intersect(
    names(ranks),
    pathway_genes
  )
  
  
  # Report ranking statistics.
  
  cat(
    "Total ranked genes:",
    length(ranks),
    "\n"
  )
  
  cat(
    "Genes overlapping Hallmark sets:",
    length(overlapping_genes),
    "\n"
  )
  
  
  # Require sufficient overlap with Hallmark gene sets.
  
  if (length(overlapping_genes) < 100) {
    
    stop(
      "Insufficient overlap between ranked genes ",
      "and Hallmark gene sets for ",
      population_name,
      ". Check gene identifier compatibility."
    )
  }
  
  
  # Report tied Wald statistics without modifying ranks.
  
  number_of_ties <- sum(
    duplicated(ranks)
  )
  
  if (number_of_ties > 0) {
    
    message(
      "Tied Wald statistics detected for ",
      population_name,
      ": ",
      number_of_ties,
      ". Original statistics will be retained."
    )
  }
  
  
  return(ranks)
}


############################################################
# 9. Run Preranked Hallmark GSEA
############################################################

run_gsea <- function(tab, population_name) {
  
  cat("\n")
  
  cat("============================================================\n")
  
  cat(
    "HALLMARK GSEA:",
    population_name,
    "\n"
  )
  
  cat("============================================================\n")
  
  
  ##########################################################
  # 9.1 Prepare Gene Ranking
  ##########################################################
  
  ranks <- make_ranks(
    tab,
    population_name
  )
  
  
  ##########################################################
  # 9.2 Save Ranked Gene Statistics
  ##########################################################
  
  ranking_df <- data.frame(
    
    gene = names(ranks),
    
    Wald_statistic = as.numeric(ranks),
    
    stringsAsFactors = FALSE
  )
  
  
  # Generate a file-safe population name.
  
  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    population_name
  )
  
  
  # Export the complete ranked gene list.
  
  write.csv(
    
    ranking_df,
    
    output_file(
      "gsea",
      paste0(
        "GSEA_Gene_Ranking_",
        safe_name,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  ##########################################################
  # 9.3 Run fgsea
  ##########################################################
  
  # Perform preranked Hallmark GSEA using Wald statistics.
  
  gsea_results <- fgsea(
    
    pathways = pathways,
    
    stats = ranks,
    
    minSize = min_size,
    
    maxSize = max_size,
    
    eps = eps_value,
    
    nproc = 1
  )
  
  
  ##########################################################
  # 9.4 Format GSEA Results
  ##########################################################
  
  # Convert fgsea output to a data frame.
  
  gsea_results <- as.data.frame(
    gsea_results
  )
  
  
  # Convert the leadingEdge list-column to character strings.
  # This allows the results to be exported as CSV.
  
  if ("leadingEdge" %in% colnames(gsea_results)) {
    
    gsea_results$leadingEdge <- vapply(
      
      gsea_results$leadingEdge,
      
      function(x) {
        paste(x, collapse = ";")
      },
      
      FUN.VALUE = character(1)
    )
  }
  
  
  # Check whether any pathways were returned.
  
  if (nrow(gsea_results) == 0) {
    
    warning(
      "No Hallmark pathways passed the size criteria for ",
      population_name
    )
    
    return(gsea_results)
  }
  
  
  # Interpret enrichment direction.
  #
  # Positive NES -> Late enrichment.
  # Negative NES -> Early enrichment.
  
  gsea_results$Stage_Enrichment <- ifelse(
    
    gsea_results$NES > 0,
    
    "Late",
    
    ifelse(
      gsea_results$NES < 0,
      "Early",
      NA
    )
  )
  
  
  # Add population information.
  
  gsea_results$Population <- population_name
  
  
  # Order pathways by adjusted P-value.
  
  gsea_results <- gsea_results[
    order(
      gsea_results$padj,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  ##########################################################
  # 9.5 Save Complete GSEA Results
  ##########################################################
  
  write.csv(
    
    gsea_results,
    
    output_file(
      "gsea",
      paste0(
        "Hallmark_GSEA_",
        safe_name,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  ##########################################################
  # 9.6 Save FDR-Significant Pathways
  ##########################################################
  
  significant_results <- gsea_results[
    
    !is.na(gsea_results$padj) &
      gsea_results$padj < fdr_cutoff,
    
    ,
    drop = FALSE
  ]
  
  
  write.csv(
    
    significant_results,
    
    output_file(
      "gsea",
      paste0(
        "Hallmark_GSEA_FDR05_",
        safe_name,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  ##########################################################
  # 9.7 Save Exploratory Nominal P-Value Results
  ##########################################################
  
  exploratory_results <- gsea_results[
    
    !is.na(gsea_results$pval) &
      gsea_results$pval < exploratory_p_cutoff,
    
    ,
    drop = FALSE
  ]
  
  
  write.csv(
    
    exploratory_results,
    
    output_file(
      "gsea",
      paste0(
        "Hallmark_GSEA_exploratory_nominalP05_",
        safe_name,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  ##########################################################
  # 9.8 Print Population-Specific Summary
  ##########################################################
  
  cat(
    "\nGSEA summary for",
    population_name,
    "\n"
  )
  
  cat(
    "Total Hallmark pathways:",
    nrow(gsea_results),
    "\n"
  )
  
  cat(
    "FDR-significant pathways:",
    nrow(significant_results),
    "\n"
  )
  
  cat(
    "Nominal P-value significant pathways:",
    nrow(exploratory_results),
    "\n"
  )
  
  
  return(gsea_results)
}


############################################################
# 10. Run GSEA for NK-like and ILC3
############################################################

# NK-like: primary population for drug repurposing.

nk_gsea <- run_gsea(
  
  tab = nk_de,
  
  population_name = "NK-like"
)


# ILC3: complementary pathway analysis.

ilc3_gsea <- run_gsea(
  
  tab = ilc3_de,
  
  population_name = "ILC3"
)


############################################################
# 11. Combine GSEA Results
############################################################

all_gsea <- rbind(
  nk_gsea,
  ilc3_gsea
)


############################################################
# 12. Save Combined GSEA Results
############################################################

write.csv(
  
  all_gsea,
  
  output_file(
    "gsea",
    "All_Populations_Hallmark_GSEA.csv"
  ),
  
  row.names = FALSE
)


############################################################
# 13. Save GSEA Results as RDS
############################################################

# Preserve population-specific result tables.

saveRDS(
  
  list(
    
    NK_like = nk_gsea,
    
    ILC3 = ilc3_gsea
  ),
  
  output_file(
    "rds",
    "12_Hallmark_GSEA_results.rds"
  )
)


############################################################
# 14. Save GSEA Analysis Parameters
############################################################

gsea_parameters <- list(
  
  dataset = "GSE184880",
  
  method = "Preranked Hallmark GSEA",
  
  ranking_metric = "DESeq2 Wald statistic",
  
  comparison = "Late vs Early",
  
  positive_NES = "Late",
  
  negative_NES = "Early",
  
  populations = c(
    "NK-like",
    "ILC3"
  ),
  
  primary_drug_repurposing_population = "NK-like",
  
  pathway_database = "MSigDB Hallmark",
  
  species = "Homo sapiens",
  
  min_size = min_size,
  
  max_size = max_size,
  
  eps = eps_value,
  
  FDR_cutoff = fdr_cutoff,
  
  exploratory_nominal_P_cutoff =
    exploratory_p_cutoff,
  
  gene_ranking_filter =
    "All genes with valid DESeq2 Wald statistics",
  
  significant_DEG_filter =
    "Not applied before GSEA"
)


# Save analysis parameters.

saveRDS(
  
  gsea_parameters,
  
  output_file(
    "gsea",
    "Hallmark_GSEA_Analysis_Parameters.rds"
  )
)


############################################################
# 15. Save Session Information
############################################################

save_session(
  "09_hallmark_gsea"
)


############################################################
# 16. Final Analysis Summary
############################################################

cat("\n")

cat("============================================================\n")

cat("SCRIPT 09: HALLMARK GSEA COMPLETED\n")

cat("============================================================\n")


# NK-like summary.

cat("\nNK-like results:\n")

cat(
  "Total pathways:",
  nrow(nk_gsea),
  "\n"
)

cat(
  "FDR-significant pathways:",
  sum(
    nk_gsea$padj < fdr_cutoff,
    na.rm = TRUE
  ),
  "\n"
)


# ILC3 summary.

cat("\nILC3 results:\n")

cat(
  "Total pathways:",
  nrow(ilc3_gsea),
  "\n"
)

cat(
  "FDR-significant pathways:",
  sum(
    ilc3_gsea$padj < fdr_cutoff,
    na.rm = TRUE
  ),
  "\n"
)


# Display output directory.

cat("\nOutput directory:\n")

cat(
  GSEA_DIR,
  "\n"
)


# Identify the primary downstream population.

cat("\nPrimary downstream population:\n")

cat(
  "NK-like -> Script 10: CMap / Drug Repurposing\n"
)


cat("\nNext step: Script 10\n")


############################################################
# END OF SCRIPT 09
############################################################
