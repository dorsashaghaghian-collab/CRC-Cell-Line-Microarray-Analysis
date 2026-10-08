# ============================================================
# CRC Cell-Line Transcriptomic Analysis
# Script 02: Expression Quality Control
# ============================================================
#
# Project:
#   GSE185055 - CRC cell-line transcriptomic analysis
#
# Purpose:
#   1. Load biological-sample metadata
#   2. Load the biological-level count matrix
#   3. Verify sample/count-matrix correspondence
#   4. Summarize library sizes
#   5. Evaluate count distributions
#   6. Apply low-count filtering
#   7. Perform variance-stabilizing transformation
#   8. Generate PCA
#   9. Generate sample-to-sample correlation analysis
#   10. Identify potential expression-level outliers
#   11. Export QC tables and diagnostic objects
#
# Experimental design:
#   4 CRC cell lines × 2 conditions × 3 biological replicates
#
# Conditions:
#   2D vs 3D
#
# Important:
#   Technical replicates were collapsed before this stage.
#   Therefore, all QC analyses are performed on 24 biological
#   samples rather than the original 192 GSM-level records.
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
  stringsAsFactors = FALSE
)

set.seed(1234)


# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "tidyverse",
  "matrixStats",
  "pheatmap"
)

missing_packages <- required_packages[
  !sapply(
    required_packages,
    requireNamespace,
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  stop(
    paste0(
      "\nThe following packages are missing:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nPlease install them before running this script."
    )
  )
}

suppressPackageStartupMessages({

  library(DESeq2)
  library(tidyverse)
  library(matrixStats)
  library(pheatmap)

})


# ------------------------------------------------------------
# 2. DIRECTORY STRUCTURE
# ------------------------------------------------------------

data_dir <- "data"
results_dir <- "results"
figures_dir <- "figures"

if (!dir.exists(data_dir)) {
  stop("Missing directory: ", data_dir)
}

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
# 3. INPUT FILES
# ------------------------------------------------------------

metadata_file <- file.path(
  data_dir,
  "sample_metadata.csv"
)

if (!file.exists(metadata_file)) {

  stop(
    "Metadata file not found:\n",
    metadata_file
  )

}


# ------------------------------------------------------------
# IMPORTANT:
# COUNT MATRIX LOCATION
# ------------------------------------------------------------
#
# The GitHub repository currently contains metadata and
# processed results, but not necessarily the full raw/processed
# count matrix.
#
# Therefore, the count matrix path must be specified here when
# running the complete analysis locally.
#
# Expected structure:
#
#   rows    = genes
#   columns = biological samples
#
# The first column should contain gene identifiers.
#
# Example:
#
# count_matrix_file <- "path/to/biological_counts.csv"
#
# ------------------------------------------------------------

count_matrix_file <- NULL


# ------------------------------------------------------------
# 4. LOAD METADATA
# ------------------------------------------------------------

message("\n============================================================")
message("Loading sample metadata")
message("============================================================")

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

required_metadata_columns <- c(
  "GSM",
  "SampleID",
  "CellLine",
  "Condition",
  "BioGroup",
  "BiologicalReplicate"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {

  stop(
    "\nMissing metadata columns:\n",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )

}

message(
  "Loaded ",
  nrow(metadata),
  " biological samples."
)


# ------------------------------------------------------------
# 5. COUNT MATRIX AVAILABILITY CHECK
# ------------------------------------------------------------

if (is.null(count_matrix_file)) {

  message("\n============================================================")
  message("COUNT MATRIX NOT CONFIGURED")
  message("============================================================")

  message(
    "\nThe metadata and QC framework are ready, ",
    "but the expression count matrix has not been linked."
  )

  message(
    "\nTo run the complete QC workflow, set:"
  )

  message(
    "\ncount_matrix_file <- ",
    '"path/to/your/count_matrix.csv"'
  )

  message(
    "\nExpected format:"
  )

  message(
    "  • rows    = genes"
  )

  message(
    "  • columns = biological samples"
  )

  message(
    "  • first column = gene identifier"
  )

  message(
    "\nThis script will stop here rather than ",
    "fabricating QC results."
  )

  quit(
    save = "no",
    status = 0
  )
}


# ------------------------------------------------------------
# 6. VERIFY COUNT MATRIX FILE
# ------------------------------------------------------------

if (!file.exists(count_matrix_file)) {

  stop(
    "\nCount matrix file does not exist:\n",
    count_matrix_file
  )

}


# ------------------------------------------------------------
# 7. LOAD COUNT MATRIX
# ------------------------------------------------------------

message("\n============================================================")
message("Loading count matrix")
message("============================================================")

counts_df <- read_csv(
  count_matrix_file,
  show_col_types = FALSE
)

message(
  "Count matrix dimensions: ",
  nrow(counts_df),
  " rows × ",
  ncol(counts_df),
  " columns."
)


# ------------------------------------------------------------
# 8. IDENTIFY GENE IDENTIFIER COLUMN
# ------------------------------------------------------------

gene_id_candidates <- c(
  "gene",
  "Gene",
  "gene_id",
  "GeneID",
  "Gene_ID",
  "ENSEMBL",
  "Ensembl",
  "GeneSymbol"
)

gene_id_column <- intersect(
  gene_id_candidates,
  colnames(counts_df)
)

if (length(gene_id_column) == 0) {

  # Assume first column is gene identifier
  gene_id_column <- colnames(counts_df)[1]

  message(
    "\nNo standard gene-ID column name detected."
  )

  message(
    "Using first column as gene identifier: ",
    gene_id_column
  )

} else {

  gene_id_column <- gene_id_column[1]

  message(
    "\nGene identifier column: ",
    gene_id_column
  )

}


# ------------------------------------------------------------
# 9. CONVERT TO MATRIX
# ------------------------------------------------------------

gene_ids <- counts_df[[gene_id_column]]

if (any(is.na(gene_ids))) {

  stop(
    "Missing gene identifiers detected in count matrix."
  )

}

if (anyDuplicated(gene_ids) > 0) {

  warning(
    "\nDuplicated gene identifiers detected."
  )

  message(
    "Duplicated identifiers: ",
    sum(
      duplicated(gene_ids)
    )
  )

}


count_data <- counts_df %>%
  select(
    -all_of(gene_id_column)
  ) %>%
  as.data.frame()

rownames(count_data) <- gene_ids


# ------------------------------------------------------------
# 10. NUMERIC COUNT VALIDATION
# ------------------------------------------------------------

non_numeric_columns <- names(
  count_data
)[
  !sapply(
    count_data,
    is.numeric
  )
]

if (length(non_numeric_columns) > 0) {

  stop(
    "\nNon-numeric expression columns detected:\n",
    paste(
      non_numeric_columns,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 11. SAMPLE MATCHING
# ------------------------------------------------------------

metadata_samples <- metadata$SampleID
count_samples <- colnames(count_data)

missing_from_counts <- setdiff(
  metadata_samples,
  count_samples
)

extra_in_counts <- setdiff(
  count_samples,
  metadata_samples
)

if (length(missing_from_counts) > 0) {

  stop(
    "\nMetadata samples missing from count matrix:\n",
    paste(
      missing_from_counts,
      collapse = ", "
    )
  )

}

if (length(extra_in_counts) > 0) {

  warning(
    "\nAdditional columns detected in count matrix:\n",
    paste(
      extra_in_counts,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 12. REORDER COUNT MATRIX
# ------------------------------------------------------------

count_data <- count_data[
  ,
  metadata_samples,
  drop = FALSE
]

if (!identical(
  colnames(count_data),
  metadata_samples
)) {

  stop(
    "\nSample order mismatch remains after reordering."
  )

}

message(
  "✓ Count matrix and metadata samples match."
)


# ------------------------------------------------------------
# 13. COUNT MATRIX BASIC STATISTICS
# ------------------------------------------------------------

message("\n============================================================")
message("Count matrix summary")
message("============================================================")

total_genes <- nrow(count_data)
total_samples <- ncol(count_data)

message(
  "Genes   : ",
  total_genes
)

message(
  "Samples : ",
  total_samples
)


# ------------------------------------------------------------
# 14. NEGATIVE / NON-INTEGER CHECK
# ------------------------------------------------------------

if (any(
  count_data < 0,
  na.rm = TRUE
)) {

  stop(
    "Negative values detected in count matrix."
  )

}

if (any(
  is.na(count_data)
)) {

  stop(
    "NA values detected in count matrix."
  )

}

non_integer_fraction <- mean(
  as.matrix(count_data) !=
    round(as.matrix(count_data))
)

if (non_integer_fraction > 0) {

  warning(
    "\nNon-integer values detected in count matrix."
  )

  message(
    "Fraction of non-integer entries: ",
    round(
      non_integer_fraction,
      5
    )
  )

}


# ------------------------------------------------------------
# 15. LIBRARY SIZE
# ------------------------------------------------------------

library_sizes <- colSums(
  count_data
)

library_size_table <- tibble(

  SampleID =
    names(library_sizes),

  Library_Size =
    as.numeric(library_sizes)

) %>%
  left_join(
    metadata,
    by = "SampleID"
  ) %>%
  arrange(
    Library_Size
  )


message("\n============================================================")
message("Library-size summary")
message("============================================================")

print(
  library_size_table
)

write_csv(
  library_size_table,
  file.path(
    results_dir,
    "library_size_summary.csv"
  )
)


# ------------------------------------------------------------
# 16. LIBRARY SIZE SUMMARY STATISTICS
# ------------------------------------------------------------

library_summary <- tibble(

  Minimum =
    min(library_sizes),

  Q1 =
    quantile(
      library_sizes,
      0.25
    ),

  Median =
    median(library_sizes),

  Mean =
    mean(library_sizes),

  Q3 =
    quantile(
      library_sizes,
      0.75
    ),

  Maximum =
    max(library_sizes)

)

print(
  library_summary
)

write_csv(
  library_summary,
  file.path(
    results_dir,
    "library_size_statistics.csv"
  )
)


# ------------------------------------------------------------
# 17. LIBRARY SIZE PLOT
# ------------------------------------------------------------

library_plot_data <- library_size_table %>%
  mutate(
    SampleID = factor(
      SampleID,
      levels = SampleID[
        order(Library_Size)
      ]
    )
  )

p_library <- ggplot(
  library_plot_data,
  aes(
    x = SampleID,
    y = Library_Size
  )
) +

  geom_col() +

  scale_y_continuous(
    labels = scales::comma
  ) +

  labs(
    title = "Library Size Across Biological Samples",
    x = "Sample",
    y = "Total Count"
  ) +

  theme_bw() +

  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
  )

ggsave(
  filename = file.path(
    figures_dir,
    "library_size.png"
  ),
  plot = p_library,
  width = 12,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 18. LOW-COUNT FILTERING
# ------------------------------------------------------------
#
# For exploratory QC and downstream DESeq2 analysis,
# genes with extremely low counts across nearly all samples
# provide little statistical information.
#
# The primary filtering rule used here is:
#
#   count >= 10 in at least 3 biological samples
#
# This threshold is deliberately documented rather than
# silently chosen.
# ------------------------------------------------------------

minimum_count <- 10
minimum_samples <- 3

keep_gene <- rowSums(
  count_data >= minimum_count
) >= minimum_samples

filtered_counts <- count_data[
  keep_gene,
  ,
  drop = FALSE
]

genes_before_filtering <- nrow(
  count_data
)

genes_after_filtering <- nrow(
  filtered_counts
)

genes_removed <- (
  genes_before_filtering -
    genes_after_filtering
)

filter_summary <- tibble(

  Genes_Before_Filtering =
    genes_before_filtering,

  Genes_After_Filtering =
    genes_after_filtering,

  Genes_Removed =
    genes_removed,

  Minimum_Count =
    minimum_count,

  Minimum_Samples =
    minimum_samples,

  Fraction_Retained =
    genes_after_filtering /
      genes_before_filtering

)

message("\n============================================================")
message("Low-count filtering")
message("============================================================")

print(
  filter_summary
)

write_csv(
  filter_summary,
  file.path(
    results_dir,
    "count_filtering_summary.csv"
  )
)


# ------------------------------------------------------------
# 19. FILTERED GENE COUNTS
# ------------------------------------------------------------

filtered_gene_counts <- tibble(

  Gene =
    rownames(filtered_counts),

  Samples_Meeting_Threshold =
    rowSums(
      filtered_counts >=
        minimum_count
    )

)

write_csv(
  filtered_gene_counts,
  file.path(
    results_dir,
    "filtered_gene_detection_summary.csv"
  )
)


# ------------------------------------------------------------
# 20. CREATE DESEQ2 OBJECT
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(

  countData =
    round(filtered_counts),

  colData =
    as.data.frame(
      metadata
    ) %>%
      column_to_rownames(
        "SampleID"
      ),

  design =
    ~ CellLine + Condition

)


# ------------------------------------------------------------
# 21. REMOVE UNUSED FACTOR LEVELS
# ------------------------------------------------------------

dds$CellLine <- droplevels(
  factor(
    dds$CellLine
  )
)

dds$Condition <- droplevels(
  factor(
    dds$Condition,
    levels = c(
      "2D",
      "3D"
    )
  )
)


# ------------------------------------------------------------
# 22. DESEQ2 SIZE FACTOR ESTIMATION
# ------------------------------------------------------------

dds <- estimateSizeFactors(
  dds
)

size_factors <- sizeFactors(
  dds
)

size_factor_table <- tibble(

  SampleID =
    names(size_factors),

  SizeFactor =
    as.numeric(size_factors)

) %>%
  left_join(
    metadata,
    by = "SampleID"
  )

write_csv(
  size_factor_table,
  file.path(
    results_dir,
    "deseq2_size_factors.csv"
  )
)


# ------------------------------------------------------------
# 23. NORMALIZED COUNTS
# ------------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

normalized_counts_export <- as.data.frame(
  normalized_counts
) %>%
  rownames_to_column(
    "Gene"
  )

# The full normalized matrix can be large.
# Export locally for reproducibility.
write_csv(
  normalized_counts_export,
  file.path(
    results_dir,
    "normalized_counts.csv"
  )
)


# ------------------------------------------------------------
# 24. VST TRANSFORMATION
# ------------------------------------------------------------

message("\n============================================================")
message("Variance-stabilizing transformation")
message("============================================================")

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vsd
)


# ------------------------------------------------------------
# 25. SAMPLE-LEVEL PCA
# ------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = c(
    "CellLine",
    "Condition"
  ),
  returnData = TRUE
)

percent_variance <- round(
  100 *
    attr(
      pca_data,
      "percentVar"
    ),
  2
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    shape = Condition
  )
) +

  geom_point(
    size = 4
  ) +

  geom_text(
    aes(
      label = name
    ),
    vjust = -0.8,
    size = 3
  ) +

  labs(
    title = "PCA of Biological Samples",
    x = paste0(
      "PC1 (",
      percent_variance[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_variance[2],
      "%)"
    )
  ) +

  facet_wrap(
    ~ CellLine
  ) +

  theme_bw()

ggsave(
  filename = file.path(
    figures_dir,
    "PCA.png"
  ),
  plot = p_pca,
  width = 11,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 26. PCA COORDINATES EXPORT
# ------------------------------------------------------------

pca_export <- pca_data %>%
  as_tibble() %>%
  rename(
    SampleID = name
  )

write_csv(
  pca_export,
  file.path(
    results_dir,
    "PCA_coordinates.csv"
  )
)


# ------------------------------------------------------------
# 27. SAMPLE-TO-SAMPLE CORRELATION
# ------------------------------------------------------------

sample_correlation <- cor(
  vst_matrix,
  method = "pearson"
)

write_csv(
  as.data.frame(
    sample_correlation
  ) %>%
    rownames_to_column(
      "SampleID"
    ),
  file.path(
    results_dir,
    "sample_correlation_matrix.csv"
  )
)


# ------------------------------------------------------------
# 28. CORRELATION HEATMAP
# ------------------------------------------------------------

annotation_col <- metadata %>%
  select(
    SampleID,
    CellLine,
    Condition
  ) %>%
  column_to_rownames(
    "SampleID"
  )

png(
  filename = file.path(
    figures_dir,
    "sample_correlation.png"
  ),
  width = 3000,
  height = 2600,
  res = 300
)

pheatmap(
  sample_correlation,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  main = "Sample-to-Sample Correlation"
)

dev.off()


# ------------------------------------------------------------
# 29. SAMPLE DISTANCE MATRIX
# ------------------------------------------------------------

sample_distance <- dist(
  t(vst_matrix)
)

sample_distance_matrix <- as.matrix(
  sample_distance
)

write_csv(
  as.data.frame(
    sample_distance_matrix
  ) %>%
    rownames_to_column(
      "SampleID"
    ),
  file.path(
    results_dir,
    "sample_distance_matrix.csv"
  )
)


# ------------------------------------------------------------
# 30. DISTANCE HEATMAP
# ------------------------------------------------------------

png(
  filename = file.path(
    figures_dir,
    "sample_distance.png"
  ),
  width = 3000,
  height = 2600,
  res = 300
)

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample-to-Sample VST Distance"
)

dev.off()


# ------------------------------------------------------------
# 31. SAMPLE-LEVEL OUTLIER SCREENING
# ------------------------------------------------------------
#
# This is a diagnostic screen, NOT an automatic deletion rule.
#
# A sample should never be removed solely because it has a
# relatively high distance or low correlation. Biological
# context and metadata must be considered.
# ------------------------------------------------------------

mean_sample_correlation <- sapply(
  seq_len(
    ncol(sample_correlation)
  ),
  function(i) {

    values <- sample_correlation[
      i,
      -i
    ]

    mean(
      values,
      na.rm = TRUE
    )
  }
)

outlier_screen <- tibble(

  SampleID =
    colnames(
      sample_correlation
    ),

  Mean_Correlation =
    mean_sample_correlation

) %>%
  left_join(
    metadata,
    by = "SampleID"
  ) %>%
  arrange(
    Mean_Correlation
  )

write_csv(
  outlier_screen,
  file.path(
    results_dir,
    "sample_outlier_screen.csv"
  )
)


# ------------------------------------------------------------
# 32. SAMPLE QC SUMMARY
# ------------------------------------------------------------

qc_summary <- library_size_table %>%

  select(
    SampleID,
    Library_Size,
    CellLine,
    Condition,
    BiologicalReplicate
  ) %>%

  left_join(
    outlier_screen %>%
      select(
        SampleID,
        Mean_Correlation
      ),
    by = "SampleID"
  ) %>%

  arrange(
    CellLine,
    Condition,
    BiologicalReplicate
  )

write_csv(
  qc_summary,
  file.path(
    results_dir,
    "sample_level_QC_summary.csv"
  )
)


# ------------------------------------------------------------
# 33. QC REPORT
# ------------------------------------------------------------

qc_report <- tibble(

  Metric = c(
    "Total genes before filtering",
    "Total genes after filtering",
    "Genes removed",
    "Total biological samples",
    "Minimum library size",
    "Median library size",
    "Maximum library size",
    "Minimum sample correlation",
    "Maximum sample correlation"
  ),

  Value = c(

    nrow(count_data),

    nrow(filtered_counts),

    nrow(count_data) -
      nrow(filtered_counts),

    ncol(count_data),

    min(library_sizes),

    median(library_sizes),

    max(library_sizes),

    min(
      sample_correlation[
        upper.tri(
          sample_correlation
        )
      ],
      na.rm = TRUE
    ),

    max(
      sample_correlation[
        upper.tri(
          sample_correlation
        )
      ],
      na.rm = TRUE
    )

  )
)

write_csv(
  qc_report,
  file.path(
    results_dir,
    "QC_summary_report.csv"
  )
)


# ------------------------------------------------------------
# 34. SAVE R OBJECTS
# ------------------------------------------------------------

saveRDS(
  dds,
  file = file.path(
    results_dir,
    "QC_DESeq2_object.rds"
  )
)

saveRDS(
  vsd,
  file = file.path(
    results_dir,
    "QC_VST_object.rds"
  )
)


# ------------------------------------------------------------
# 35. SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    results_dir,
    "QC_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 36. FINAL MESSAGE
# ------------------------------------------------------------

message("\n============================================================")
message("QUALITY CONTROL COMPLETED")
message("============================================================")

message(
  "\nGenes before filtering : ",
  nrow(count_data)
)

message(
  "Genes after filtering  : ",
  nrow(filtered_counts)
)

message(
  "Biological samples     : ",
  ncol(count_data)
)

message(
  "Minimum library size   : ",
  min(library_sizes)
)

message(
  "Median library size    : ",
  median(library_sizes)
)

message(
  "Maximum library size   : ",
  max(library_sizes)
)

message("\nGenerated QC outputs include:")
message("  ✓ library_size_summary.csv")
message("  ✓ count_filtering_summary.csv")
message("  ✓ normalized_counts.csv")
message("  ✓ PCA_coordinates.csv")
message("  ✓ sample_correlation_matrix.csv")
message("  ✓ sample_distance_matrix.csv")
message("  ✓ sample_outlier_screen.csv")
message("  ✓ sample_level_QC_summary.csv")
message("  ✓ QC_summary_report.csv")
message("  ✓ PCA.png")
message("  ✓ sample_correlation.png")
message("  ✓ sample_distance.png")

message("\nIMPORTANT:")
message(
  "Outlier screening is diagnostic only; ",
  "samples are NOT automatically removed."
)

message("============================================================\n")
