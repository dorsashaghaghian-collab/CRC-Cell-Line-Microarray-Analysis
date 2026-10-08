# ============================================================
# CRC Cell-Line Transcriptomic Analysis
# Script 04: Chemokine–Receptor Panel Analysis
# ============================================================
#
# Project:
#   GSE185055 - CRC cell-line transcriptomic analysis
#
# Objective:
#   Characterize the expression response of a predefined
#   chemokine/receptor panel between 3D and 2D culture
#   conditions.
#
# Biological panel:
#
#   CXCL9  ─┐
#   CXCL10 ├── CXCR3
#   CXCL11 ─┘
#
#   CCL3 ─┐
#   CCL4 ├── CCR5
#   CCL5 ─┘
#
#   CCL2 ─┐
#   CCL7 ├── CCR2
#   CCL8 ─┘
#
#   CXCL12 ── CXCR4
#   CXCL16 ── CXCR6
#   CX3CL1 ── CX3CR1
#
# Analysis levels:
#   1. Global 3D vs 2D
#   2. Within-cell-line 3D vs 2D
#   3. Panel-level biological interpretation
#
# Important:
#   NA values in differential-expression results indicate that
#   a statistical result was not available for that gene.
#   They must NOT automatically be interpreted as biological
#   absence of expression.
#
# Author:
#   Dorsa Shaghaghian
#
# ============================================================


# ------------------------------------------------------------
# 0. CLEAN SESSION
# ------------------------------------------------------------

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

set.seed(1234)


# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "readr",
  "dplyr",
  "tibble",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  stop(
    paste0(
      "\nMissing required packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nInstall the missing packages before running ",
      "this script."
    )
  )

}

suppressPackageStartupMessages({

  library(tidyverse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)

})


# ------------------------------------------------------------
# 2. PROJECT DIRECTORIES
# ------------------------------------------------------------

data_dir <- "data"
results_dir <- "results"
figures_dir <- "figures"

if (!dir.exists(results_dir)) {

  dir.create(
    results_dir,
    recursive = TRUE
  )

}

if (!dir.exists(figures_dir)) {

  dir.create(
    figures_dir,
    recursive = TRUE
  )

}


# ------------------------------------------------------------
# 3. ANALYSIS PARAMETERS
# ------------------------------------------------------------

alpha <- 0.05

lfc_threshold <- 1

panel_name <- "CRC_chemokine_receptor_panel"


# ------------------------------------------------------------
# 4. DEFINE BIOLOGICAL PANEL
# ------------------------------------------------------------
#
# The panel is defined explicitly rather than being inferred
# from significant results. This prevents selection bias.
# ------------------------------------------------------------

chemokine_panel <- tribble(

  ~Gene,     ~Role,       ~Axis,

  "CXCL9",   "Ligand",    "CXCR3",
  "CXCL10",  "Ligand",    "CXCR3",
  "CXCL11",  "Ligand",    "CXCR3",

  "CCL3",    "Ligand",    "CCR5",
  "CCL4",    "Ligand",    "CCR5",
  "CCL5",    "Ligand",    "CCR5",

  "CCL2",    "Ligand",    "CCR2",
  "CCL7",    "Ligand",    "CCR2",
  "CCL8",    "Ligand",    "CCR2",

  "CXCL12",  "Ligand",    "CXCR4",
  "CXCR4",   "Receptor",  "CXCR4",

  "CXCL16",  "Ligand",    "CXCR6",
  "CXCR6",   "Receptor",  "CXCR6",

  "CX3CL1",  "Ligand",    "CX3CR1",
  "CX3CR1",  "Receptor",  "CX3CR1",

  "CXCR3",   "Receptor",  "CXCR3",
  "CCR5",    "Receptor",  "CCR5",
  "CCR2",    "Receptor",  "CCR2"

)


# ------------------------------------------------------------
# 5. VALIDATE PANEL
# ------------------------------------------------------------

panel_genes <- unique(
  chemokine_panel$Gene
)

expected_panel_size <- 18

if (
  length(panel_genes) !=
  expected_panel_size
) {

  stop(
    "\nUnexpected number of unique genes in panel.\n",
    "Expected: ",
    expected_panel_size,
    "\nObserved: ",
    length(panel_genes)
  )

}

if (
  anyDuplicated(
    chemokine_panel$Gene
  ) > 0
) {

  stop(
    "\nDuplicate genes detected in panel definition."
  )

}


# ------------------------------------------------------------
# 6. SAVE PANEL DEFINITION
# ------------------------------------------------------------

write_csv(

  chemokine_panel,

  file.path(
    results_dir,
    "chemokine_receptor_panel_definition.csv"
  )

)


# ------------------------------------------------------------
# 7. DEFINE LIGAND–RECEPTOR AXES
# ------------------------------------------------------------

axis_definition <- tribble(

  ~Axis,      ~Ligands,                     ~Receptor,

  "CXCR3",
  "CXCL9;CXCL10;CXCL11",
  "CXCR3",

  "CCR5",
  "CCL3;CCL4;CCL5",
  "CCR5",

  "CCR2",
  "CCL2;CCL7;CCL8",
  "CCR2",

  "CXCR4",
  "CXCL12",
  "CXCR4",

  "CXCR6",
  "CXCL16",
  "CXCR6",

  "CX3CR1",
  "CX3CL1",
  "CX3CR1"

)


write_csv(

  axis_definition,

  file.path(
    results_dir,
    "chemokine_receptor_axes.csv"
  )

)


# ------------------------------------------------------------
# 8. INPUT RESULT FILES
# ------------------------------------------------------------

global_file <- file.path(
  results_dir,
  "global_3D_vs_2D.csv"
)

within_file <- file.path(
  results_dir,
  "within_cell_line_3D_vs_2D.csv"
)


if (!file.exists(global_file)) {

  stop(
    "\nGlobal DE result file not found:\n",
    global_file,
    "\n\nRun Script 03 first."
  )

}

if (!file.exists(within_file)) {

  stop(
    "\nWithin-cell-line DE result file not found:\n",
    within_file,
    "\n\nRun Script 03 first."
  )

}


# ------------------------------------------------------------
# 9. LOAD DIFFERENTIAL-EXPRESSION RESULTS
# ------------------------------------------------------------

message("\n============================================================")
message("Loading differential-expression results")
message("============================================================")

global_de <- read_csv(
  global_file,
  show_col_types = FALSE
)

within_de <- read_csv(
  within_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 10. VALIDATE REQUIRED DE COLUMNS
# ------------------------------------------------------------

required_de_columns <- c(
  "Gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

missing_global_columns <- setdiff(
  required_de_columns,
  colnames(global_de)
)

if (length(missing_global_columns) > 0) {

  stop(
    "\nGlobal DE results are missing required columns:\n",
    paste(
      missing_global_columns,
      collapse = ", "
    )
  )

}

missing_within_columns <- setdiff(
  required_de_columns,
  colnames(within_de)
)

if (length(missing_within_columns) > 0) {

  stop(
    "\nWithin-cell-line DE results are missing required columns:\n",
    paste(
      missing_within_columns,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 11. CHECK PANEL COVERAGE
# ------------------------------------------------------------

global_panel_genes <- intersect(
  panel_genes,
  global_de$Gene
)

within_panel_genes <- intersect(
  panel_genes,
  within_de$Gene
)

global_missing_genes <- setdiff(
  panel_genes,
  global_de$Gene
)

within_missing_genes <- setdiff(
  panel_genes,
  within_de$Gene
)


# ------------------------------------------------------------
# 12. PANEL COVERAGE SUMMARY
# ------------------------------------------------------------

coverage_summary <- tibble(

  Analysis = c(
    "Global",
    "Within-cell-line"
  ),

  Panel_Genes =
    c(
      length(global_panel_genes),
      length(within_panel_genes)
    ),

  Missing_From_Result =
    c(
      length(global_missing_genes),
      length(within_missing_genes)
    ),

  Total_Panel_Genes =
    length(panel_genes)

)

print(
  coverage_summary
)

write_csv(

  coverage_summary,

  file.path(
    results_dir,
    "chemokine_panel_coverage_summary.csv"
  )

)


# ------------------------------------------------------------
# 13. GLOBAL PANEL RESULTS
# ------------------------------------------------------------

global_panel <- global_de %>%

  filter(
    Gene %in% panel_genes
  ) %>%

  left_join(
    chemokine_panel,
    by = "Gene"
  ) %>%

  mutate(

    FDR_Significant =
      !is.na(padj) &
      padj < alpha,

    Effect_Size_Significant =
      !is.na(log2FoldChange) &
      abs(log2FoldChange) >=
        lfc_threshold,

    Significant =
      FDR_Significant &
      Effect_Size_Significant,

    Direction = case_when(

      Significant &
        log2FoldChange > 0 ~
        "Up_in_3D",

      Significant &
        log2FoldChange < 0 ~
        "Down_in_3D",

      TRUE ~
        "Not_significant_or_unavailable"

    )

  ) %>%

  select(

    Gene,
    Role,
    Axis,

    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,

    FDR_Significant,
    Effect_Size_Significant,
    Significant,
    Direction

  ) %>%

  arrange(
    Axis,
    desc(
      Role == "Receptor"
    ),
    Gene
  )


# ------------------------------------------------------------
# 14. PRESERVE PANEL ORDER
# ------------------------------------------------------------

global_panel <- global_panel %>%

  mutate(

    Gene = factor(
      Gene,
      levels = panel_genes
    )

  ) %>%

  arrange(
    Gene
  ) %>%

  mutate(
    Gene = as.character(Gene)
  )


# ------------------------------------------------------------
# 15. SAVE GLOBAL PANEL
# ------------------------------------------------------------

write_csv(

  global_panel,

  file.path(
    results_dir,
    "chemokine_receptor_panel.csv"
  )

)


# ------------------------------------------------------------
# 16. GLOBAL PANEL SUMMARY
# ------------------------------------------------------------

global_panel_summary <- global_panel %>%

  summarise(

    Panel_Genes =
      n(),

    Genes_With_DE_Result =
      sum(
        !is.na(padj)
      ),

    FDR_Significant =
      sum(
        FDR_Significant,
        na.rm = TRUE
      ),

    Significant_Both =
      sum(
        Significant,
        na.rm = TRUE
      ),

    Up_in_3D =
      sum(
        Direction == "Up_in_3D"
      ),

    Down_in_3D =
      sum(
        Direction == "Down_in_3D"
      )

  )

print(
  global_panel_summary
)

write_csv(

  global_panel_summary,

  file.path(
    results_dir,
    "global_chemokine_panel_summary.csv"
  )

)


# ------------------------------------------------------------
# 17. WITHIN-CELL-LINE PANEL RESULTS
# ------------------------------------------------------------

within_panel <- within_de %>%

  filter(
    Gene %in% panel_genes
  ) %>%

  left_join(
    chemokine_panel,
    by = "Gene"
  ) %>%

  mutate(

    FDR_Significant =
      !is.na(padj) &
      padj < alpha,

    Effect_Size_Significant =
      !is.na(log2FoldChange) &
      abs(log2FoldChange) >=
        lfc_threshold,

    Significant =
      FDR_Significant &
      Effect_Size_Significant,

    Direction = case_when(

      Significant &
        log2FoldChange > 0 ~
        "Up_in_3D",

      Significant &
        log2FoldChange < 0 ~
        "Down_in_3D",

      TRUE ~
        "Not_significant_or_unavailable"

    )

  ) %>%

  select(

    CellLine,

    Gene,
    Role,
    Axis,

    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,

    FDR_Significant,
    Effect_Size_Significant,
    Significant,
    Direction

  ) %>%

  arrange(
    CellLine,
    Axis,
    Gene
  )


# ------------------------------------------------------------
# 18. SAVE WITHIN-CELL-LINE PANEL
# ------------------------------------------------------------

write_csv(

  within_panel,

  file.path(
    results_dir,
    "within_cell_line_chemokine_panel.csv"
  )

)


# ------------------------------------------------------------
# 19. CELL-LINE PANEL SUMMARY
# ------------------------------------------------------------

within_panel_summary <- within_panel %>%

  group_by(
    CellLine
  ) %>%

  summarise(

    Panel_Genes =
      n(),

    Genes_With_DE_Result =
      sum(
        !is.na(padj)
      ),

    FDR_Significant =
      sum(
        FDR_Significant,
        na.rm = TRUE
      ),

    Significant_Both =
      sum(
        Significant,
        na.rm = TRUE
      ),

    Up_in_3D =
      sum(
        Direction == "Up_in_3D"
      ),

    Down_in_3D =
      sum(
        Direction == "Down_in_3D"
      ),

    .groups = "drop"

  )

print(
  within_panel_summary
)

write_csv(

  within_panel_summary,

  file.path(
    results_dir,
    "within_cell_line_chemokine_panel_summary.csv"
  )

)


# ------------------------------------------------------------
# 20. AXIS-LEVEL SUMMARY
# ------------------------------------------------------------
#
# This summarizes the predefined ligand/receptor axes.
#
# IMPORTANT:
# This is a descriptive summary and does not imply that all
# components of an axis are simultaneously active.
# ------------------------------------------------------------

axis_summary <- within_panel %>%

  group_by(
    CellLine,
    Axis
  ) %>%

  summarise(

    Number_of_Panel_Genes =
      n(),

    Number_with_DE_result =
      sum(
        !is.na(padj)
      ),

    Significant_Genes =
      sum(
        Significant,
        na.rm = TRUE
      ),

    Significant_Up =
      sum(
        Direction == "Up_in_3D"
      ),

    Significant_Down =
      sum(
        Direction == "Down_in_3D"
      ),

    .groups = "drop"

  )

write_csv(

  axis_summary,

  file.path(
    results_dir,
    "chemokine_receptor_axis_summary.csv"
  )

)


# ------------------------------------------------------------
# 21. GLOBAL AXIS SUMMARY
# ------------------------------------------------------------

global_axis_summary <- global_panel %>%

  group_by(
    Axis
  ) %>%

  summarise(

    Number_of_Panel_Genes =
      n(),

    Genes_with_DE_result =
      sum(
        !is.na(padj)
      ),

    Significant_Genes =
      sum(
        Significant,
        na.rm = TRUE
      ),

    Significant_Up =
      sum(
        Direction == "Up_in_3D"
      ),

    Significant_Down =
      sum(
        Direction == "Down_in_3D"
      ),

    .groups = "drop"

  )

write_csv(

  global_axis_summary,

  file.path(
    results_dir,
    "global_chemokine_receptor_axis_summary.csv"
  )

)


# ------------------------------------------------------------
# 22. PANEL DIRECTION MATRIX
# ------------------------------------------------------------

direction_matrix <- within_panel %>%

  select(
    CellLine,
    Gene,
    Direction
  ) %>%

  pivot_wider(

    names_from =
      CellLine,

    values_from =
      Direction

  )

write_csv(

  direction_matrix,

  file.path(
    results_dir,
    "chemokine_panel_direction_matrix.csv"
  )

)


# ------------------------------------------------------------
# 23. PANEL LOG2FC MATRIX
# ------------------------------------------------------------

log2fc_matrix <- within_panel %>%

  select(
    CellLine,
    Gene,
    log2FoldChange
  ) %>%

  pivot_wider(

    names_from =
      CellLine,

    values_from =
      log2FoldChange

  )

write_csv(

  log2fc_matrix,

  file.path(
    results_dir,
    "chemokine_panel_log2FC_matrix.csv"
  )

)


# ------------------------------------------------------------
# 24. PANEL FDR MATRIX
# ------------------------------------------------------------

fdr_matrix <- within_panel %>%

  select(
    CellLine,
    Gene,
    padj
  ) %>%

  pivot_wider(

    names_from =
      CellLine,

    values_from =
      padj

  )

write_csv(

  fdr_matrix,

  file.path(
    results_dir,
    "chemokine_panel_FDR_matrix.csv"
  )

)


# ------------------------------------------------------------
# 25. GLOBAL PANEL VOLCANO-LIKE DATA
# ------------------------------------------------------------
#
# The panel is not used to define genome-wide significance.
# Instead, the genome-wide DE results are retained and panel
# genes are highlighted as a separate annotation layer.
# ------------------------------------------------------------

global_panel_plot_data <- global_de %>%

  mutate(

    Panel_Gene =
      Gene %in% panel_genes,

    NegLog10_Padj =
      -log10(
        pmax(
          padj,
          .Machine$double.xmin
        )
      )

  )


# ------------------------------------------------------------
# 26. GLOBAL PANEL HIGHLIGHT PLOT
# ------------------------------------------------------------

p_panel_volcano <- ggplot(

  global_panel_plot_data,

  aes(
    x = log2FoldChange,
    y = NegLog10_Padj
  )

) +

  geom_point(
    alpha = 0.25,
    size = 1
  ) +

  geom_point(
    data =
      subset(
        global_panel_plot_data,
        Panel_Gene
      ),
    size = 2.8
  ) +

  geom_vline(
    xintercept =
      c(
        -lfc_threshold,
        lfc_threshold
      ),
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept =
      -log10(alpha),
    linetype = "dashed"
  ) +

  geom_text(
    data =
      subset(
        global_panel_plot_data,
        Panel_Gene &
        !is.na(log2FoldChange) &
        !is.na(NegLog10_Padj)
      ),
    aes(
      label = Gene
    ),
    vjust = -0.7,
    size = 3,
    check_overlap = TRUE
  ) +

  labs(

    title =
      "Genome-wide Differential Expression with Chemokine/Receptor Panel",

    subtitle =
      "Panel genes highlighted within the global 3D vs 2D analysis",

    x =
      "log2 Fold Change (3D vs 2D)",

    y =
      "-log10 adjusted P value"

  ) +

  theme_bw()


ggsave(

  filename =
    file.path(
      figures_dir,
      "chemokine_panel_highlight_volcano.png"
    ),

  plot =
    p_panel_volcano,

  width =
    11,

  height =
    8,

  dpi =
    300

)


# ------------------------------------------------------------
# 27. PANEL LOG2FC HEATMAP DATA
# ------------------------------------------------------------

heatmap_data <- within_panel %>%

  select(
    Gene,
    CellLine,
    log2FoldChange
  ) %>%

  pivot_wider(

    names_from =
      CellLine,

    values_from =
      log2FoldChange

  ) %>%

  arrange(
    match(
      Gene,
      panel_genes
    )
  )


write_csv(

  heatmap_data,

  file.path(
    results_dir,
    "chemokine_panel_heatmap_data.csv"
  )

)


# ------------------------------------------------------------
# 28. PANEL HEATMAP
# ------------------------------------------------------------

heatmap_matrix <- heatmap_data %>%

  column_to_rownames(
    "Gene"
  ) %>%

  as.matrix()


# Keep the matrix numeric
storage.mode(
  heatmap_matrix
) <- "numeric"


# ------------------------------------------------------------
# 29. SAVE HEATMAP IMAGE
# ------------------------------------------------------------

if (
  nrow(heatmap_matrix) > 0 &&
  ncol(heatmap_matrix) > 0
) {

  png(

    filename =
      file.path(
        figures_dir,
        "chemokine_receptor_heatmap.png"
      ),

    width =
      2600,

    height =
      3000,

    res =
      300

  )

  pheatmap::pheatmap(

    heatmap_matrix,

    cluster_rows = FALSE,

    cluster_cols = FALSE,

    na_col = "grey90",

    border_color = NA,

    main =
      "Chemokine/Receptor Panel: log2FC (3D vs 2D)"

  )

  dev.off()

}


# ------------------------------------------------------------
# 30. IDENTIFY KEY PANEL FINDINGS
# ------------------------------------------------------------

key_global_findings <- global_panel %>%

  filter(
    Significant
  ) %>%

  arrange(
    padj
  )

write_csv(

  key_global_findings,

  file.path(
    results_dir,
    "key_global_chemokine_receptor_findings.csv"
  )

)


# ------------------------------------------------------------
# 31. IDENTIFY CELL-LINE-SPECIFIC FINDINGS
# ------------------------------------------------------------

key_cellline_findings <- within_panel %>%

  filter(
    Significant
  ) %>%

  arrange(
    CellLine,
    padj
  )

write_csv(

  key_cellline_findings,

  file.path(
    results_dir,
    "key_cell_line_chemokine_receptor_findings.csv"
  )

)


# ------------------------------------------------------------
# 32. CREATE A HUMAN-READABLE SUMMARY
# ------------------------------------------------------------

summary_lines <- c(

  "CRC CHEMOKINE/RECEPTOR PANEL ANALYSIS",
  "======================================",
  "",
  paste0(
    "Panel genes: ",
    length(panel_genes)
  ),

  paste0(
    "Global panel genes with DE results: ",
    sum(
      !is.na(
        global_panel$padj
      )
    )
  ),

  paste0(
    "Global panel genes significant at FDR < ",
    alpha,
    " and |log2FC| >= ",
    lfc_threshold,
    ": ",
    sum(
      global_panel$Significant,
      na.rm = TRUE
    )
  ),

  "",

  "Axes:",
  paste(
    axis_definition$Axis,
    collapse = ", "
  ),

  "",

  "Interpretation note:",
  "NA differential-expression statistics indicate that a",
  "statistical result was unavailable and should not be",
  "interpreted as biological absence of expression."

)

writeLines(

  summary_lines,

  con =
    file.path(
      results_dir,
      "chemokine_receptor_analysis_summary.txt"
    )

)


# ------------------------------------------------------------
# 33. SAVE SESSION INFORMATION
# ------------------------------------------------------------

capture.output(

  sessionInfo(),

  file =
    file.path(
      results_dir,
      "chemokine_receptor_sessionInfo.txt"
    )

)


# ------------------------------------------------------------
# 34. FINAL VALIDATION
# ------------------------------------------------------------

message("\n============================================================")
message("CHEMOKINE/RECEPTOR PANEL ANALYSIS COMPLETED")
message("============================================================")

message(
  "\nPanel genes defined : ",
  length(panel_genes)
)

message(
  "Global DE results   : ",
  length(global_panel_genes),
  " panel genes found"
)

message(
  "Global significant  : ",
  sum(
    global_panel$Significant,
    na.rm = TRUE
  )
)

message(
  "\nCell-line analyses:"
)

for (
  cell_line in unique(
    within_panel$CellLine
  )
) {

  n_sig <- within_panel %>%

    filter(
      CellLine == cell_line
    ) %>%

    summarise(
      n = sum(
        Significant,
        na.rm = TRUE
      )
    ) %>%

    pull(n)

  message(
    "  ✓ ",
    cell_line,
    ": ",
    n_sig,
    " significant panel genes"
  )

}

message(
  "\nSignificance criterion:"
)

message(
  "  FDR < ",
  alpha,
  " AND |log2FC| >= ",
  lfc_threshold
)

message(
  "\nImportant:"
)

message(
  "NA values are treated as unavailable statistical results, ",
  "not as evidence of biological absence."
)

message("============================================================\n")
