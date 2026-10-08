# ============================================================
# 05_visualization.R
# CRC Cell-Line Transcriptomic Analysis
#
# Purpose:
#   Generate publication-quality visualizations for:
#   1. PCA
#   2. Sample correlation
#   3. Global 3D vs 2D volcano plot
#   4. Chemokine/receptor panel heatmap
#
# Study:
#   GSE185055
#
# Design:
#   4 CRC cell lines × 2 conditions × 3 biological replicates
#
# Author:
#   Dorsa Shaghaghian
#
# Reproducibility:
#   All figures are generated from validated analysis outputs.
# ============================================================


# ------------------------------------------------------------
# 0. CLEAN SESSION AND OPTIONS
# ------------------------------------------------------------

rm(list = ls())

options(
  stringsAsFactors = FALSE,
  scipen = 999
)


# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "pheatmap",
  "RColorBrewer"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following required packages are not installed:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nPlease install them before running this script."
    )
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})


# ------------------------------------------------------------
# 2. DIRECTORY STRUCTURE
# ------------------------------------------------------------

project_dir <- "E:/CRC_GEO_EXPRESSION_DATA/CRC_CELL_LINE_GITHUB_PACKAGE"

results_dir <- file.path(
  project_dir,
  "results"
)

figures_dir <- file.path(
  project_dir,
  "figures"
)

data_dir <- file.path(
  project_dir,
  "data"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------

check_file <- function(path, description) {

  if (!file.exists(path)) {
    stop(
      paste0(
        "\nRequired file not found: ",
        description,
        "\nExpected path:\n",
        path,
        "\n\nPlease run the preceding analysis script first."
      )
    )
  }

  invisible(path)
}


check_columns <- function(
    data,
    required_columns,
    dataset_name
) {

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "\nMissing required columns in ",
        dataset_name,
        ":\n",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}


save_figure <- function(
    filename,
    width,
    height
) {

  output_path <- file.path(
    figures_dir,
    filename
  )

  ggsave(
    filename = output_path,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  message(
    "Saved figure: ",
    output_path
  )
}


# ------------------------------------------------------------
# 4. INPUT FILES
# ------------------------------------------------------------

metadata_file <- file.path(
  data_dir,
  "sample_metadata.csv"
)

global_de_file <- file.path(
  results_dir,
  "global_3D_vs_2D.csv"
)

within_de_file <- file.path(
  results_dir,
  "within_cell_line_3D_vs_2D.csv"
)

chemokine_panel_file <- file.path(
  results_dir,
  "chemokine_receptor_panel.csv"
)

pca_coordinates_file <- file.path(
  results_dir,
  "PCA_coordinates.csv"
)

correlation_matrix_file <- file.path(
  results_dir,
  "sample_correlation_matrix.csv"
)


# ------------------------------------------------------------
# 5. VERIFY REQUIRED INPUTS
# ------------------------------------------------------------

check_file(
  metadata_file,
  "sample metadata"
)

check_file(
  global_de_file,
  "global differential-expression results"
)

check_file(
  within_de_file,
  "within-cell-line differential-expression results"
)

check_file(
  chemokine_panel_file,
  "chemokine/receptor panel results"
)


# ------------------------------------------------------------
# 6. LOAD DATA
# ------------------------------------------------------------

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

global_de <- read_csv(
  global_de_file,
  show_col_types = FALSE
)

within_de <- read_csv(
  within_de_file,
  show_col_types = FALSE
)

chemokine_panel <- read_csv(
  chemokine_panel_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 7. VALIDATE METADATA
# ------------------------------------------------------------

check_columns(
  metadata,
  c(
    "GSM",
    "SampleID",
    "CellLine",
    "Condition",
    "BioGroup",
    "BiologicalReplicate"
  ),
  "sample_metadata.csv"
)

if (nrow(metadata) != 24) {

  stop(
    paste0(
      "Expected 24 final biological samples, but found ",
      nrow(metadata),
      "."
    )
  )
}

if (anyDuplicated(metadata$SampleID) > 0) {

  stop(
    "Duplicate SampleID values detected in metadata."
  )
}

expected_cell_lines <- c(
  "HCT116",
  "HT29",
  "LS174T",
  "LS513"
)

expected_conditions <- c(
  "2D",
  "3D"
)

if (!all(expected_cell_lines %in% metadata$CellLine)) {

  stop(
    "Metadata does not contain all expected CRC cell lines."
  )
}

if (!all(expected_conditions %in% metadata$Condition)) {

  stop(
    "Metadata does not contain both 2D and 3D conditions."
  )
}


# ------------------------------------------------------------
# 8. VALIDATE GLOBAL DE RESULTS
# ------------------------------------------------------------

check_columns(
  global_de,
  c(
    "gene",
    "log2FoldChange",
    "padj"
  ),
  "global_3D_vs_2D.csv"
)


# ------------------------------------------------------------
# 9. VALIDATE WITHIN-CELL-LINE RESULTS
# ------------------------------------------------------------

check_columns(
  within_de,
  c(
    "CellLine",
    "gene",
    "log2FoldChange",
    "padj"
  ),
  "within_cell_line_3D_vs_2D.csv"
)


# ------------------------------------------------------------
# 10. VALIDATE CHEMOKINE/RECEPTOR RESULTS
# ------------------------------------------------------------

check_columns(
  chemokine_panel,
  c(
    "Gene",
    "log2FoldChange",
    "padj"
  ),
  "chemokine_receptor_panel.csv"
)


# ------------------------------------------------------------
# 11. GLOBAL VISUALIZATION PARAMETERS
# ------------------------------------------------------------

alpha_threshold <- 0.05

log2fc_threshold <- 1

panel_genes <- c(
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "CCL3",
  "CCL4",
  "CCL5",
  "CCL2",
  "CCL7",
  "CCL8",
  "CXCL12",
  "CXCR4",
  "CXCL16",
  "CXCR6",
  "CX3CL1",
  "CX3CR1",
  "CXCR3",
  "CCR5",
  "CCR2"
)


# ============================================================
# FIGURE 1
# PCA
# ============================================================

if (file.exists(pca_coordinates_file)) {

  pca_data <- read_csv(
    pca_coordinates_file,
    show_col_types = FALSE
  )

  check_columns(
    pca_data,
    c(
      "SampleID",
      "PC1",
      "PC2"
    ),
    "PCA_coordinates.csv"
  )

  pca_data <- pca_data %>%
    left_join(
      metadata %>%
        select(
          SampleID,
          CellLine,
          Condition,
          BioGroup
        ),
      by = "SampleID"
    )

  if (any(is.na(pca_data$CellLine))) {

    stop(
      "Some PCA samples could not be matched to metadata."
    )
  }

  pca_plot <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      shape = CellLine,
      fill = Condition
    )
  ) +

    geom_point(
      size = 4,
      stroke = 0.8,
      colour = "black"
    ) +

    theme_classic(
      base_size = 13
    ) +

    labs(
      title = "PCA of CRC Cell-Line Transcriptomes",
      subtitle = "GSE185055",
      x = "Principal Component 1",
      y = "Principal Component 2",
      shape = "Cell line",
      fill = "Condition"
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(
        size = 11
      ),
      legend.position = "right"
    )

  save_figure(
    "PCA.png",
    width = 9,
    height = 6.5
  )

  print(pca_plot)

} else {

  message(
    "\nPCA_coordinates.csv was not found."
  )

  message(
    "PCA.png will not be regenerated by this script."
  )

  message(
    "Run 02_quality_control.R first if PCA coordinates are required."
  )
}


# ============================================================
# FIGURE 2
# SAMPLE CORRELATION
# ============================================================

if (file.exists(correlation_matrix_file)) {

  correlation_matrix <- read_csv(
    correlation_matrix_file,
    show_col_types = FALSE
  )

  correlation_matrix <- as.data.frame(
    correlation_matrix
  )

  if (ncol(correlation_matrix) < 2) {

    stop(
      "Sample correlation matrix does not contain enough columns."
    )
  }

  rownames(correlation_matrix) <-
    correlation_matrix[[1]]

  correlation_matrix[[1]] <- NULL

  correlation_matrix <- as.matrix(
    correlation_matrix
  )

  if (!is.numeric(correlation_matrix)) {

    correlation_matrix <- apply(
      correlation_matrix,
      2,
      as.numeric
    )

    rownames(correlation_matrix) <-
      rownames(
        as.data.frame(
          read_csv(
            correlation_matrix_file,
            show_col_types = FALSE
          )
        )
      )
  }

  annotation_samples <- data.frame(
    CellLine = metadata$CellLine,
    Condition = metadata$Condition
  )

  rownames(annotation_samples) <-
    metadata$SampleID

  common_samples <- intersect(
    colnames(correlation_matrix),
    rownames(annotation_samples)
  )

  if (length(common_samples) >= 2) {

    correlation_matrix <- correlation_matrix[
      common_samples,
      common_samples,
      drop = FALSE
    ]

    annotation_samples <- annotation_samples[
      common_samples,
      ,
      drop = FALSE
    ]

    png(
      filename = file.path(
        figures_dir,
        "sample_correlation.png"
      ),
      width = 2600,
      height = 2400,
      res = 300
    )

    pheatmap(
      correlation_matrix,
      annotation_col = annotation_samples,
      annotation_row = annotation_samples,
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      clustering_method = "complete",
      border_color = NA,
      main = "Sample-to-Sample Correlation"
    )

    dev.off()

    message(
      "Saved figure: ",
      file.path(
        figures_dir,
        "sample_correlation.png"
      )
    )

  } else {

    message(
      "Insufficient matched samples for correlation heatmap."
    )
  }

} else {

  message(
    "\nSample correlation matrix was not found."
  )

  message(
    "sample_correlation.png will not be regenerated."
  )
}


# ============================================================
# FIGURE 3
# GLOBAL 3D VS 2D VOLCANO PLOT
# ============================================================

volcano_data <- global_de %>%

  mutate(
    neg_log10_padj = if_else(
      !is.na(padj) & padj > 0,
      -log10(padj),
      NA_real_
    ),

    Significance = case_when(

      !is.na(padj) &
        padj < alpha_threshold &
        log2FoldChange >= log2fc_threshold ~
        "Upregulated in 3D",

      !is.na(padj) &
        padj < alpha_threshold &
        log2FoldChange <= -log2fc_threshold ~
        "Downregulated in 3D",

      TRUE ~
        "Not significant"
    ),

    Panel = gene %in% panel_genes
  )


volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj
  )
) +

  geom_point(
    aes(
      alpha = Significance
    ),
    size = 1.5
  ) +

  geom_vline(
    xintercept = c(
      -log2fc_threshold,
      log2fc_threshold
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = -log10(alpha_threshold),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_point(
    data = volcano_data %>%
      filter(
        Panel,
        !is.na(neg_log10_padj),
        !is.na(log2FoldChange)
      ),
    size = 2.8,
    shape = 21,
    stroke = 0.8
  ) +

  geom_text(
    data = volcano_data %>%
      filter(
        Panel,
        !is.na(neg_log10_padj),
        !is.na(log2FoldChange)
      ),
    aes(
      label = gene
    ),
    size = 3,
    vjust = -0.8,
    check_overlap = TRUE
  ) +

  scale_alpha_manual(
    values = c(
      "Upregulated in 3D" = 0.8,
      "Downregulated in 3D" = 0.8,
      "Not significant" = 0.25
    )
  ) +

  theme_classic(
    base_size = 13
  ) +

  labs(
    title = "Global Differential Expression",
    subtitle = "3D vs 2D CRC cell-line models",
    x = "log2 Fold Change",
    y = expression(-log[10]("adjusted P-value")),
    alpha = "Classification"
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    legend.position = "right"
  )


ggsave(
  filename = file.path(
    figures_dir,
    "volcano_3D_vs_2D.png"
  ),
  plot = volcano_plot,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "Saved figure: ",
  file.path(
    figures_dir,
    "volcano_3D_vs_2D.png"
  )
)

print(volcano_plot)


# ============================================================
# FIGURE 4
# CHEMOKINE/RECEPTOR PANEL HEATMAP
# ============================================================

panel_heatmap_data <- chemokine_panel %>%

  filter(
    Gene %in% panel_genes
  ) %>%

  select(
    Gene,
    CellLine,
    log2FoldChange
  ) %>%

  distinct()


if ("CellLine" %in% colnames(panel_heatmap_data)) {

  panel_heatmap_matrix <- panel_heatmap_data %>%

    pivot_wider(
      names_from = CellLine,
      values_from = log2FoldChange
    ) %>%

    column_to_rownames(
      "Gene"
    ) %>%

    as.matrix()


  missing_panel_genes <- setdiff(
    panel_genes,
    rownames(panel_heatmap_matrix)
  )


  if (length(missing_panel_genes) > 0) {

    missing_rows <- matrix(
      NA_real_,
      nrow = length(missing_panel_genes),
      ncol = ncol(panel_heatmap_matrix),
      dimnames = list(
        missing_panel_genes,
        colnames(panel_heatmap_matrix)
      )
    )

    panel_heatmap_matrix <- rbind(
      panel_heatmap_matrix,
      missing_rows
    )
  }


  panel_heatmap_matrix <- panel_heatmap_matrix[
    panel_genes,
    ,
    drop = FALSE
  ]


  png(
    filename = file.path(
      figures_dir,
      "chemokine_receptor_heatmap.png"
    ),
    width = 2800,
    height = 3600,
    res = 300
  )

  pheatmap(
    panel_heatmap_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    na_col = "grey90",
    main = "Chemokine/Receptor Panel\n3D vs 2D log2 Fold Change",
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 11
  )

  dev.off()

  message(
    "Saved figure: ",
    file.path(
      figures_dir,
      "chemokine_receptor_heatmap.png"
    )
  )

} else {

  stop(
    paste0(
      "CellLine column is required in ",
      "chemokine_receptor_panel.csv ",
      "to generate the panel heatmap."
    )
  )
}


# ============================================================
# 12. FIGURE MANIFEST
# ============================================================

figure_manifest <- tibble(

  Figure = c(
    "PCA.png",
    "sample_correlation.png",
    "volcano_3D_vs_2D.png",
    "chemokine_receptor_heatmap.png"
  ),

  Description = c(
    "Principal component analysis of biological samples",
    "Sample-to-sample expression correlation",
    "Global differential expression: 3D versus 2D",
    "Chemokine/receptor panel across CRC cell lines"
  ),

  Main_Input = c(
    "PCA_coordinates.csv + sample_metadata.csv",
    "sample_correlation_matrix.csv + sample_metadata.csv",
    "global_3D_vs_2D.csv",
    "chemokine_receptor_panel.csv"
  ),

  Resolution = rep(
    "300 dpi",
    4
  )
)


write_csv(
  figure_manifest,
  file.path(
    results_dir,
    "figure_manifest.csv"
  )
)


# ============================================================
# 13. SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    results_dir,
    "visualization_sessionInfo.txt"
  )
)


# ============================================================
# 14. FINAL VALIDATION
# ============================================================

expected_figures <- c(
  "volcano_3D_vs_2D.png",
  "chemokine_receptor_heatmap.png"
)

existing_figures <- file.exists(
  file.path(
    figures_dir,
    expected_figures
  )
)

if (!all(existing_figures)) {

  stop(
    paste0(
      "\nVisualization validation failed.\n",
      "Missing required figure(s):\n",
      paste(
        expected_figures[!existing_figures],
        collapse = "\n"
      )
    )
  )
}


message("\n==============================================")
message("VISUALIZATION PIPELINE COMPLETED SUCCESSFULLY")
message("==============================================")
message("Figures directory:")
message(figures_dir)
message("")
message("Generated/validated outputs:")
message("- PCA.png")
message("- sample_correlation.png")
message("- volcano_3D_vs_2D.png")
message("- chemokine_receptor_heatmap.png")
message("- figure_manifest.csv")
message("- visualization_sessionInfo.txt")
message("==============================================")
