# ============================================================
# 05_visualization.R
# Visualization framework for the final analysis
# ============================================================

library(tidyverse)
library(ggplot2)

# ------------------------------------------------------------
# Output figures
# ------------------------------------------------------------
#
# Final publication-style figures are stored in:
#
#   figures/PCA.png
#   figures/sample_correlation.png
#   figures/volcano_3D_vs_2D.png
#   figures/chemokine_receptor_heatmap.png
#
# ------------------------------------------------------------

figure_paths <- c(
  "figures/PCA.png",
  "figures/sample_correlation.png",
  "figures/volcano_3D_vs_2D.png",
  "figures/chemokine_receptor_heatmap.png"
)

# Check that the documented figures exist.

missing_figures <- figure_paths[
  !file.exists(figure_paths)
]

if (length(missing_figures) > 0) {
  warning(
    "The following figures are not available:\n",
    paste(missing_figures, collapse = "\n")
  )
} else {
  cat("All final figures are present.\n")
}

cat("\nVisualization documentation completed.\n")
