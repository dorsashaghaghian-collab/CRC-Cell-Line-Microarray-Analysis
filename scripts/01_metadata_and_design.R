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
    paste(colnames(metadata), collapse = ", ")
  )
}

message("All required metadata columns are present.")


# ------------------------------------------------------------
# 5. BASIC DATA-INTEGRITY CHECKS
# ------------------------------------------------------------

message("\n------------------------------------------------------------")
message("Checking metadata integrity")
message("------------------------------------------------------------")


# Number of rows
if (nrow(metadata) == 0) {
  stop("Metadata contains zero rows.")
}


# Check missing values in required columns
missing_summary <- metadata %>%
  summarise(
    across(
      all_of(required_columns),
      ~ sum(is.na(.))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Missing_Count"
  ) %>%
  mutate(
    Missing = Missing_Count > 0
  )

print(missing_summary)

write_csv(
  missing_summary,
  file.path(results_dir, "metadata_missing_values.csv")
)

if (any(missing_summary$Missing)) {
  warning(
    "\nMissing values detected in one or more required metadata fields."
  )
} else {
  message("No missing values detected in required metadata fields.")
}


# ------------------------------------------------------------
# 6. DUPLICATE SAMPLE CHECK
# ------------------------------------------------------------

message("\n------------------------------------------------------------")
message("Checking for duplicate biological samples")
message("------------------------------------------------------------")


duplicate_sample_ids <- metadata %>%
  count(SampleID, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_sample_ids) > 0) {

  print(duplicate_sample_ids)

  stop(
    "\nDuplicate SampleID values detected.\n",
    "Each biological sample must have a unique SampleID."
  )

} else {

  message("No duplicate SampleID values detected.")
}


# ------------------------------------------------------------
# 7. DUPLICATE GSM CHECK
# ------------------------------------------------------------

duplicate_gsm <- metadata %>%
  count(GSM, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_gsm) > 0) {

  print(duplicate_gsm)

  stop(
    "\nDuplicate GSM identifiers detected."
  )

} else {

  message("No duplicate GSM identifiers detected.")
}


# ------------------------------------------------------------
# 8. VALIDATE EXPECTED CELL LINES
# ------------------------------------------------------------

expected_cell_lines <- c(
  "HCT116",
  "HT29",
  "LS174T",
  "LS513"
)

observed_cell_lines <- sort(
  unique(metadata$CellLine)
)

message("\nObserved cell lines:")
print(observed_cell_lines)

unexpected_cell_lines <- setdiff(
  observed_cell_lines,
  expected_cell_lines
)

missing_expected_cell_lines <- setdiff(
  expected_cell_lines,
  observed_cell_lines
)

if (length(unexpected_cell_lines) > 0) {

  warning(
    "\nUnexpected cell lines detected:\n",
    paste(unexpected_cell_lines, collapse = ", ")
  )
}

if (length(missing_expected_cell_lines) > 0) {

  warning(
    "\nExpected cell lines missing from metadata:\n",
    paste(missing_expected_cell_lines, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 9. VALIDATE EXPERIMENTAL CONDITIONS
# ------------------------------------------------------------

expected_conditions <- c(
  "2D",
  "3D"
)

observed_conditions <- sort(
  unique(metadata$Condition)
)

message("\nObserved experimental conditions:")
print(observed_conditions)

unexpected_conditions <- setdiff(
  observed_conditions,
  expected_conditions
)

if (length(unexpected_conditions) > 0) {

  stop(
    "\nUnexpected experimental condition(s) detected:\n",
    paste(unexpected_conditions, collapse = ", "),
    "\n\nExpected conditions: 2D and 3D."
  )
}


# ------------------------------------------------------------
# 10. STANDARDIZE FACTOR LEVELS
# ------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    CellLine = factor(
      CellLine,
      levels = expected_cell_lines
    ),

    Condition = factor(
      Condition,
      levels = expected_conditions
    )
  )


# ------------------------------------------------------------
# 11. CREATE ANALYSIS GROUP
# ------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    Group = paste(
      CellLine,
      Condition,
      sep = "_"
    )
  )

message("\nAnalysis groups:")
print(
  metadata %>%
    count(Group)
)


# ------------------------------------------------------------
# 12. CELL LINE × CONDITION DESIGN TABLE
# ------------------------------------------------------------

design_table <- metadata %>%
  count(
    CellLine,
    Condition,
    name = "Biological_Sample_Count"
  ) %>%
  arrange(
    CellLine,
    Condition
  )

message("\n============================================================")
message("Experimental design")
message("============================================================")

print(design_table)

write_csv(
  design_table,
  file.path(
    results_dir,
    "design_cellline_condition_counts.csv"
  )
)


# ------------------------------------------------------------
# 13. EXPECTED REPLICATE VALIDATION
# ------------------------------------------------------------

expected_replicates <- 3

replicate_validation <- metadata %>%
  group_by(
    CellLine,
    Condition
  ) %>%
  summarise(
    Biological_Replicates =
      n_distinct(BiologicalReplicate),

    Samples =
      n(),

    .groups = "drop"
  ) %>%
  mutate(
    Expected_Replicates = expected_replicates,

    Replicate_Count_Correct =
      Biological_Replicates == expected_replicates,

    Sample_Count_Correct =
      Samples == expected_replicates,

    Design_Valid =
      Replicate_Count_Correct &
      Sample_Count_Correct
  )

message("\nBiological replicate validation:")
print(replicate_validation)

write_csv(
  replicate_validation,
  file.path(
    results_dir,
    "biological_replicate_validation.csv"
  )
)

if (!all(replicate_validation$Design_Valid)) {

  stop(
    "\nExperimental design validation failed.\n",
    "Every cell-line × condition combination should contain ",
    expected_replicates,
    " biological replicates."
  )
}

message(
  "\n✓ All cell-line × condition groups contain exactly ",
  expected_replicates,
  " biological replicates."
)


# ------------------------------------------------------------
# 14. BIOLOGICAL REPLICATE IDENTIFIER CHECK
# ------------------------------------------------------------

replicate_values <- metadata %>%
  distinct(
    BiologicalReplicate
  ) %>%
  arrange(
    BiologicalReplicate
  )

message("\nObserved biological replicate identifiers:")
print(replicate_values)


# Check whether replicate numbers are exactly 1, 2, 3
expected_replicate_ids <- c(1, 2, 3)

observed_replicate_ids <- sort(
  unique(
    metadata$BiologicalReplicate
  )
)

if (!identical(
  as.numeric(observed_replicate_ids),
  expected_replicate_ids
)) {

  warning(
    "\nBiological replicate identifiers are not exactly 1, 2, 3."
  )
}


# ------------------------------------------------------------
# 15. CHECK CELL LINE × CONDITION × REPLICATE UNIQUENESS
# ------------------------------------------------------------

replicate_duplicates <- metadata %>%
  count(
    CellLine,
    Condition,
    BiologicalReplicate,
    name = "n"
  ) %>%
  filter(n > 1)

if (nrow(replicate_duplicates) > 0) {

  print(replicate_duplicates)

  stop(
    "\nDuplicate biological replicate assignments detected."
  )

} else {

  message(
    "✓ Each biological replicate is uniquely assigned."
  )
}


# ------------------------------------------------------------
# 16. BIOGROUP VALIDATION
# ------------------------------------------------------------

if ("BioGroup" %in% colnames(metadata)) {

  expected_biogroup <- paste(
    metadata$CellLine,
    metadata$Condition,
    sep = "_"
  )

  biogroup_matches <- metadata$BioGroup == expected_biogroup

  if (all(biogroup_matches)) {

    message(
      "✓ BioGroup values are consistent with CellLine × Condition."
    )

  } else {

    warning(
      "\nSome BioGroup values do not match ",
      "CellLine × Condition."
    )

    inconsistent_biogroups <- metadata %>%
      mutate(
        Expected_BioGroup =
          paste(CellLine, Condition, sep = "_")
      ) %>%
      filter(
        BioGroup != Expected_BioGroup
      ) %>%
      select(
        GSM,
        SampleID,
        CellLine,
        Condition,
        BioGroup,
        Expected_BioGroup
      )

    print(inconsistent_biogroups)

    write_csv(
      inconsistent_biogroups,
      file.path(
        results_dir,
        "inconsistent_biogroup_assignments.csv"
      )
    )
  }
}


# ------------------------------------------------------------
# 17. SAMPLE-LEVEL DESIGN TABLE
# ------------------------------------------------------------

sample_design <- metadata %>%
  select(
    GSM,
    SampleID,
    CellLine,
    Condition,
    Group,
    BioGroup,
    BiologicalReplicate
  ) %>%
  arrange(
    CellLine,
    Condition,
    BiologicalReplicate
  )

write_csv(
  sample_design,
  file.path(
    results_dir,
    "analysis_ready_sample_design.csv"
  )
)


# ------------------------------------------------------------
# 18. SAMPLE DISTRIBUTION SUMMARY
# ------------------------------------------------------------

sample_summary <- metadata %>%
  summarise(
    Total_Biological_Samples =
      n(),

    Number_of_Cell_Lines =
      n_distinct(CellLine),

    Number_of_Conditions =
      n_distinct(Condition),

    Number_of_Biological_Replicates =
      n_distinct(BiologicalReplicate),

    Number_of_Experimental_Groups =
      n_distinct(Group)
  )

message("\n============================================================")
message("Study summary")
message("============================================================")

print(sample_summary)

write_csv(
  sample_summary,
  file.path(
    results_dir,
    "study_design_summary.csv"
  )
)


# ------------------------------------------------------------
# 19. EXPECTED STUDY STRUCTURE CHECK
# ------------------------------------------------------------

expected_total_samples <- 24
expected_groups <- 8

actual_total_samples <- nrow(metadata)
actual_groups <- n_distinct(metadata$Group)

if (actual_total_samples != expected_total_samples) {

  stop(
    "\nUnexpected number of biological samples.\n",
    "Expected: ",
    expected_total_samples,
    "\nObserved: ",
    actual_total_samples
  )
}

if (actual_groups != expected_groups) {

  stop(
    "\nUnexpected number of experimental groups.\n",
    "Expected: ",
    expected_groups,
    "\nObserved: ",
    actual_groups
  )
}

message(
  "\n✓ Study contains the expected ",
  expected_total_samples,
  " biological samples across ",
  expected_groups,
  " experimental groups."
)


# ------------------------------------------------------------
# 20. GLOBAL CONDITION DISTRIBUTION
# ------------------------------------------------------------

condition_summary <- metadata %>%
  count(
    Condition,
    name = "Biological_Samples"
  )

message("\nSamples by condition:")
print(condition_summary)

write_csv(
  condition_summary,
  file.path(
    results_dir,
    "condition_summary.csv"
  )
)


# ------------------------------------------------------------
# 21. CELL-LINE DISTRIBUTION
# ------------------------------------------------------------

cellline_summary <- metadata %>%
  count(
    CellLine,
    name = "Biological_Samples"
  )

message("\nSamples by cell line:")
print(cellline_summary)

write_csv(
  cellline_summary,
  file.path(
    results_dir,
    "cell_line_summary.csv"
  )
)


# ------------------------------------------------------------
# 22. FINAL METADATA OBJECT
# ------------------------------------------------------------

# Convert factors to explicit character columns before export.
# This makes the CSV portable across R, Python, Excel and other
# analysis environments.

metadata_export <- metadata %>%
  mutate(
    CellLine = as.character(CellLine),
    Condition = as.character(Condition)
  )

write_csv(
  metadata_export,
  file.path(
    results_dir,
    "validated_metadata.csv"
  )
)


# ------------------------------------------------------------
# 23. SESSION INFORMATION
# ------------------------------------------------------------

session_file <- file.path(
  results_dir,
  "metadata_validation_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_file
)


# ------------------------------------------------------------
# 24. FINAL VALIDATION REPORT
# ------------------------------------------------------------

validation_status <- tibble(
  Check = c(
    "Metadata file exists",
    "Required columns present",
    "No missing required metadata",
    "Unique SampleID",
    "Unique GSM",
    "Expected cell lines",
    "Expected conditions",
    "Three biological replicates per group",
    "Unique biological replicate assignments",
    "BioGroup consistency",
    "24 biological samples",
    "8 experimental groups"
  ),

  Status = c(
    file.exists(metadata_file),
    length(missing_columns) == 0,
    !any(missing_summary$Missing),
    nrow(duplicate_sample_ids) == 0,
    nrow(duplicate_gsm) == 0,
    length(unexpected_cell_lines) == 0,
    length(unexpected_conditions) == 0,
    all(replicate_validation$Design_Valid),
    nrow(replicate_duplicates) == 0,
    all(biogroup_matches),
    actual_total_samples == expected_total_samples,
    actual_groups == expected_groups
  )
)

write_csv(
  validation_status,
  file.path(
    results_dir,
    "metadata_validation_report.csv"
  )
)


# ------------------------------------------------------------
# 25. FINAL STATUS
# ------------------------------------------------------------

message("\n")
message("============================================================")
message("METADATA VALIDATION COMPLETED SUCCESSFULLY")
message("============================================================")

message(
  "Biological samples : ",
  actual_total_samples
)

message(
  "Cell lines         : ",
  n_distinct(metadata$CellLine)
)

message(
  "Conditions         : ",
  n_distinct(metadata$Condition)
)

message(
  "Experimental groups: ",
  actual_groups
)

message(
  "Replicates/group   : ",
  expected_replicates
)

message("\nGenerated validation files:")

generated_files <- c(
  "metadata_missing_values.csv",
  "design_cellline_condition_counts.csv",
  "biological_replicate_validation.csv",
  "analysis_ready_sample_design.csv",
  "study_design_summary.csv",
  "condition_summary.csv",
  "cell_line_summary.csv",
  "validated_metadata.csv",
  "metadata_validation_sessionInfo.txt",
  "metadata_validation_report.csv"
)

for (file in generated_files) {
  message(
    "  ✓ ",
    file.path(results_dir, file)
  )
}

message("\nNo downstream analysis should proceed if this validation fails.")
message("============================================================\n")
