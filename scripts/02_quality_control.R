# ============================================================
# 02_quality_control.R
# Quality control and exploratory sample assessment
# ============================================================

library(tidyverse)
library(ggplot2)

# ------------------------------------------------------------
# NOTE
# ------------------------------------------------------------
# The raw GEO expression matrix is intentionally not stored
# in this repository.
#
# This script documents the QC workflow used for the final
# biological-level dataset.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 1. Load sample metadata
# ------------------------------------------------------------

metadata <- read.csv(
  "data/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Inspect experimental structure
# ------------------------------------------------------------

print(
  metadata %>%
    count(CellLine, Condition)
)

# ------------------------------------------------------------
# 3. QC outputs
# ------------------------------------------------------------
#
# Final QC figures generated from the processed expression
# matrix are provided in:
#
#   figures/PCA.png
#   figures/sample_correlation.png
#
# The corresponding analysis evaluates:
#
#   - library-size consistency
#   - sample-level structure
#   - replicate coherence
#   - potential outliers
#
# ------------------------------------------------------------

cat("\nQC documentation completed.\n")
cat("See figures/PCA.png and figures/sample_correlation.png\n")
