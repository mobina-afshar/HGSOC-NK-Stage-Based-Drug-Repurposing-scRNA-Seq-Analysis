
############################################################
# GSE184880 - HGSOC scRNA-seq Analysis
# Script 01: Project Setup
#
# Purpose:
# Initialize the project directory structure, load the shared
# configuration, validate sample metadata, and prepare output
# directories for downstream analysis.
#
# Input:
# - config.yaml
#
# Output:
# - results/metadata/sample_metadata.csv
# - Project directory structure
# - logs/01_setup_sessionInfo.txt
#
# IMPORTANT:
# - Run this script from the repository root.
# - Existing input files and analysis results are not deleted.
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

# Define the project root and main directories.

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

# Create the main project directories.

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

# Set the project-wide random seed.

set.seed(cfg$project$seed)


############################################################
# 2. Define Shared Helper Functions
############################################################

# Generate output file paths and create parent directories
# automatically when needed.

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


# Check whether all required R packages are installed.

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
# 3. Create Project Directory Structure
############################################################

# Define the directories required throughout the analysis.

required_dirs <- c(
  
  # Input and reference data
  RAW_DIR,
  file.path(
    PROJECT_ROOT,
    cfg$paths$reference
  ),
  
  # Quality control and preprocessing
  file.path(RESULTS_DIR, "qc"),
  file.path(RESULTS_DIR, "preprocessing"),
  
  # Annotation and cell population analysis
  file.path(RESULTS_DIR, "annotation"),
  file.path(RESULTS_DIR, "tnk_ilc"),
  
  # Pseudobulk and differential expression
  file.path(RESULTS_DIR, "pseudobulk"),
  file.path(RESULTS_DIR, "differential_expression"),
  
  # Pathway analysis and drug repurposing
  file.path(RESULTS_DIR, "gsea"),
  file.path(RESULTS_DIR, "cmap"),
  file.path(RESULTS_DIR, "drug_annotation"),
  
  # R objects and log files
  RDS_DIR,
  LOGS_DIR
)

# Create all required directories.

invisible(
  lapply(
    required_dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)


############################################################
# 4. Validate Sample Metadata
############################################################

# Extract sample information from config.yaml.

sample_metadata <- data.frame(
  
  sample = unlist(cfg$samples$sample),
  
  GSM = unlist(cfg$samples$gsm),
  
  stage = unlist(cfg$samples$stage),
  
  stage_group = unlist(cfg$samples$stage_group),
  
  stringsAsFactors = FALSE
)


# Check sample and GSM identifier uniqueness.

if (
  anyDuplicated(sample_metadata$sample) ||
  anyDuplicated(sample_metadata$GSM)
) {
  
  stop(
    "Sample identifiers in config.yaml must be unique."
  )
}


# Validate the stage group labels.

if (
  !all(
    sample_metadata$stage_group %in% c("Early", "Late")
  )
) {
  
  stop(
    "stage_group must contain only 'Early' and 'Late'."
  )
}


# Save the validated sample metadata.

write.csv(
  sample_metadata,
  out_dir("metadata", "sample_metadata.csv"),
  row.names = FALSE
)


############################################################
# 5. Report Project Initialization
############################################################

cat("\n")
cat("============================================================\n")
cat("PROJECT INITIALIZATION COMPLETED\n")
cat("============================================================\n")

cat("Project root:", PROJECT_ROOT, "\n")

cat("Raw 10x input directory:", RAW_DIR, "\n")

cat(
  "Number of samples:",
  nrow(sample_metadata),
  "\n"
)

cat(
  "Expected input files per sample:\n",
  "<GSM>_<CancerN>.matrix.mtx.gz\n",
  "<GSM>_<CancerN>.barcodes.tsv.gz\n",
  "<GSM>_<CancerN>.genes.tsv.gz\n",
  sep = ""
)

cat(
  "Sample metadata saved to: ",
  file.path(RESULTS_DIR, "metadata", "sample_metadata.csv"),
  "\n",
  sep = ""
)

cat(
  "Existing results and input files have not been deleted.\n"
)


############################################################
# 6. Save Session Information
############################################################

save_session("01_setup")


############################################################
# END OF SCRIPT 01
############################################################
