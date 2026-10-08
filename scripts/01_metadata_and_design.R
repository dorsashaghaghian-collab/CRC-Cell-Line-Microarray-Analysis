# ============================================================
# 01_metadata_and_design.R
# CRC Cell-Line Transcriptomic Analysis
#
# Purpose:
#   Define and validate the biological-level experimental design.
#
# Dataset:
#   GSE185055
#
# Design:
#   4 colorectal cancer cell lines
#   2D vs 3D growth
#   3 biological replicates per cell-line × condition
# ============================================================

# Required packages
library(tidyverse)

# ------------------------------------------------------------
# 1. Load curated metadata
# ------------------------------------------------------------

metadata <- read.csv(
  "data/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Basic structure checks
# ------------------------------------------------------------

cat("Number of biological samples:", nrow(metadata), "\n")
cat("Number of cell lines:", n_distinct(metadata$CellLine), "\n")
cat("Conditions:", paste(unique(metadata$Condition), collapse = ", "), "\n")

# ------------------------------------------------------------
# 3. Validate experimental design
# ------------------------------------------------------------

design_table <- metadata %>%
  count(CellLine, Condition, name = "BiologicalReplicates")

print(design_table)

# Expected design:
#
# CellLine × Condition
# --------------------
# HCT116   × 2D
# HCT116   × 3D
# HT29     × 2D
# HT29     × 3D
# LS174T   × 2D
# LS174T   × 3D
# LS513    × 2D
# LS513    × 3D
#
# Each combination should contain 3 biological replicates.

stopifnot(all(design_table$BiologicalReplicates == 3))

# ------------------------------------------------------------
# 4. Create reproducible sample grouping
# ------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    Group = interaction(
      CellLine,
      Condition,
      sep = "_"
    )
  )

# ------------------------------------------------------------
# 5. Save validated design table
# ------------------------------------------------------------

write.csv(
  design_table,
  "results/design_validation.csv",
  row.names = FALSE
)

cat("\nDesign validation completed successfully.\n")
