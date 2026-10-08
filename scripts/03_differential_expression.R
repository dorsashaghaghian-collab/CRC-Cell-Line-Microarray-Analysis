# ============================================================
# 03_differential_expression.R
# Differential expression framework
#
# Primary comparison:
#   3D vs 2D
#
# Analysis levels:
#   1. Global comparison
#   2. Cell-line-specific comparisons
# ============================================================

library(DESeq2)
library(tidyverse)

# ------------------------------------------------------------
# 1. Load biological-level metadata
# ------------------------------------------------------------

metadata <- read.csv(
  "data/sample_metadata.csv",
  stringsAsFactors = FALSE
)

metadata$Condition <- factor(
  metadata$Condition,
  levels = c("2D", "3D")
)

# ------------------------------------------------------------
# 2. Expression matrix
# ------------------------------------------------------------
#
# The raw count matrix is not stored in this repository.
# It should be supplied as a matrix with:
#
#   rows    = genes
#   columns = biological samples
#
# Column names must match metadata$SampleID.
# ------------------------------------------------------------

# counts <- read.csv(
#   "path/to/count_matrix.csv",
#   row.names = 1,
#   check.names = FALSE
# )

# ------------------------------------------------------------
# 3. Construct DESeq2 object
# ------------------------------------------------------------

# dds <- DESeqDataSetFromMatrix(
#   countData = counts,
#   colData = metadata,
#   design = ~ Condition
# )

# ------------------------------------------------------------
# 4. Run differential expression
# ------------------------------------------------------------

# dds <- DESeq(dds)

# Global 3D vs 2D:
#
# res_global <- results(
#   dds,
#   contrast = c("Condition", "3D", "2D")
# )

# ------------------------------------------------------------
# 5. Significance criteria
# ------------------------------------------------------------

# sig_global <- as.data.frame(res_global) %>%
#   rownames_to_column("Gene") %>%
#   filter(
#     !is.na(padj),
#     padj < 0.05,
#     abs(log2FoldChange) >= 1
#   )

# ------------------------------------------------------------
# 6. Cell-line-specific analyses
# ------------------------------------------------------------
#
# For each cell line, the same 3D vs 2D comparison is performed
# using only samples belonging to that cell line.
#
# The final repository contains the resulting summary table:
#
#   results/within_cell_line_3D_vs_2D.csv
# ------------------------------------------------------------

cat("\nDifferential-expression workflow documented.\n")
cat("Primary contrast: 3D vs 2D\n")
cat("Threshold: FDR < 0.05 and |log2FC| >= 1\n")
