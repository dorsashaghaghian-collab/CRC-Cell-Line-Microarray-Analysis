# ============================================================
# CRC Cell-Line Transcriptomic Analysis
# Script 01: Metadata and Experimental Design Validation
# ============================================================
#
# Project:
#   GSE185055 - CRC cell-line transcriptomic analysis
#
# Purpose:
#   1. Load and validate the biological-sample metadata
#   2. Confirm the experimental design
#   3. Check biological replicate structure
#   4. Detect missing values and duplicate samples
#   5. Generate analysis-ready design variables
#   6. Export validation tables for reproducibility
#
# Experimental design:
#   Cell lines : HCT116, HT29, LS174T, LS513
#   Conditions : 2D, 3D
#   Biological replicates : 3 per cell-line × condition
#
# Important:
#   The analysis is performed at the biological-sample level.
#   Technical replicates were collapsed before downstream
#   differential-expression analysis.
#
# Author:
#   Dorsa Shaghaghian
#
# ============================================================


# ------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------

# Clean workspace
rm(list = ls())

# Avoid automatic conversion of strings
options(stringsAsFactors = FALSE)

# Reproducibility
set.seed(1234)


# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "readr",
  "dplyr",
  "tibble"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "\nThe following required packages are missing:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them before running this script."
    )
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(tibble)
})


# ------------------------------------------------------------
# 2. PROJECT DIRECTORIES
# ------------------------------------------------------------

# The script is designed to be executed from the repository root.
#
# Expected repository structure:
#
# CRC-Cell-Line-Microarray-Analysis/
# ├── data/
# ├── docs/
# ├── figures/
# ├── results/
# ├── scripts/
# └── README.md

data_dir <- "data"
results_dir <- "results"

metadata_file <- file.path(
  data_dir,
  "sample_metadata.csv"
)

if (!dir.exists(data_dir)) {
  stop("Directory not found: ", data_dir)
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found:\n",
    metadata_file,
    "\n\nRun this script from the repository root."
  )
}


# ------------------------------------------------------------
# 3. LOAD METADATA
# ------------------------------------------------------------

message("\n============================================================")
message("Loading biological-sample metadata")
message("============================================================")

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

message(
  "Metadata loaded successfully: ",
  nrow(metadata),
  " rows × ",
  ncol(metadata),
  " columns."
)


# ------------------------------------------------------------
# 4. REQUIRED COLUMN VALIDATION
# ------------------------------------------------------------

required_columns <- c(
  "GSM",
  "SampleID",
  "CellLine",
  "Condition",
  "BioGroup",
  "BiologicalReplicate"
)

missing_columns <- setdiff(
  required_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "\nMissing required metadata columns:\n",
    paste(missing_columns, collapse = ", "),
    "\n\nAvailable columns:\n",
    paste
