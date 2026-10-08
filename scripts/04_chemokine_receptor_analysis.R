# ============================================================
# 04_chemokine_receptor_analysis.R
# Chemokine–receptor focused analysis
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Define the predefined gene panel
# ------------------------------------------------------------

chemokine_receptor_panel <- c(
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CXCR3",
  "CCL3",
  "CCL4",
  "CCL5",
  "CCR5",
  "CCL2",
  "CCL7",
  "CCL8",
  "CCR2",
  "CXCL12",
  "CXCR4",
  "CXCL16",
  "CXCR6",
  "CX3CL1",
  "CX3CR1"
)

# ------------------------------------------------------------
# 2. Define ligand–receptor axes
# ------------------------------------------------------------

chemokine_axes <- tibble(
  Axis = c(
    "CXCR3",
    "CCR5",
    "CCR2",
    "CXCR4",
    "CXCR6",
    "CX3CR1"
  ),
  Ligands = c(
    "CXCL9, CXCL10, CXCL11",
    "CCL3, CCL4, CCL5",
    "CCL2, CCL7, CCL8",
    "CXCL12",
    "CXCL16",
    "CX3CL1"
  )
)

# ------------------------------------------------------------
# 3. Load final panel results
# ------------------------------------------------------------

panel_results <- read.csv(
  "results/chemokine_receptor_panel.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 4. Report panel coverage
# ------------------------------------------------------------

cat(
  "Number of predefined panel genes:",
  length(chemokine_receptor_panel),
  "\n"
)

cat(
  "Number of result rows:",
  nrow(panel_results),
  "\n"
)

# ------------------------------------------------------------
# 5. Save panel definition
# ------------------------------------------------------------

write.csv(
  chemokine_axes,
  "results/chemokine_receptor_axes.csv",
  row.names = FALSE
)

cat("\nChemokine–receptor panel analysis documented.\n")
