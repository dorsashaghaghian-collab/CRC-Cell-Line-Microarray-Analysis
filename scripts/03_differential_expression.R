# ============================================================
# CRC Cell-Line Transcriptomic Analysis
# Script 03: Differential Expression Analysis
# ============================================================
#
# Project:
#   GSE185055 - CRC cell-line transcriptomic analysis
#
# Objective:
#   Identify genes differentially expressed between 3D and 2D
#   culture conditions at:
#
#   1. Global level across all four CRC cell lines
#   2. Individual cell-line level
#
# Experimental design:
#   4 cell lines × 2 conditions × 3 biological replicates
#
# Cell lines:
#   HCT116, HT29, LS174T, LS513
#
# Conditions:
#   2D, 3D
#
# Statistical framework:
#   DESeq2
#
# Significance criteria:
#   Adjusted P value (FDR) < 0.05
#   AND
#   absolute log2 fold-change >= 1
#
# Important:
#   Technical replicates were collapsed before differential
#   expression analysis. Statistical inference is therefore
#   performed at the biological-replicate level.
#
# Author:
#   Dorsa Shaghaghian
#
# ============================================================


# ------------------------------------------------------------
# 0. CLEAN SESSION AND GLOBAL OPTIONS
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
  "DESeq2",
  "tidyverse",
  "readr",
  "dplyr",
  "tibble"
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
      "\nMissing required R packages:\n",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nPlease install the missing packages before ",
      "running this analysis."
    )
  )
}

suppressPackageStartupMessages({

  library(DESeq2)
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(tibble)

})


# ------------------------------------------------------------
# 2. PROJECT DIRECTORIES
# ------------------------------------------------------------

data_dir <- "data"
results_dir <- "results"

if (!dir.exists(data_dir)) {

  stop(
    "\nData directory not found:\n",
    data_dir,
    "\n\nRun this script from the repository root."
  )

}

if (!dir.exists(results_dir)) {

  dir.create(
    results_dir,
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
    "\nMetadata file not found:\n",
    metadata_file
  )

}


# ------------------------------------------------------------
# COUNT MATRIX CONFIGURATION
# ------------------------------------------------------------
#
# The complete count matrix is intentionally not assumed to
# exist inside the public GitHub repository.
#
# When running the full analysis locally, replace NULL with
# the path to the final biological-level count matrix.
#
# Expected structure:
#
#   rows    = genes
#   columns = biological samples
#   first column = gene identifier
#
# Example:
#
# count_matrix_file <- "E:/CRC_GEO_EXPRESSION_DATA/....csv"
#
# ------------------------------------------------------------

count_matrix_file <- NULL


# ------------------------------------------------------------
# 4. ANALYSIS PARAMETERS
# ------------------------------------------------------------

# Reference condition
reference_condition <- "2D"

# Experimental condition
experimental_condition <- "3D"

# Multiple-testing threshold
alpha <- 0.05

# Minimum biologically relevant effect size
lfc_threshold <- 1

# Independent filtering is handled by DESeq2.
# No arbitrary post-hoc filtering is applied to adjusted P values.


# ------------------------------------------------------------
# 5. LOAD METADATA
# ------------------------------------------------------------

message("\n============================================================")
message("Loading metadata")
message("============================================================")

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 6. VALIDATE METADATA STRUCTURE
# ------------------------------------------------------------

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
    "\nRequired metadata columns are missing:\n",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 7. VALIDATE SAMPLE IDENTIFIERS
# ------------------------------------------------------------

if (anyDuplicated(metadata$SampleID) > 0) {

  duplicated_samples <- metadata %>%
    count(
      SampleID,
      name = "n"
    ) %>%
    filter(
      n > 1
    )

  print(duplicated_samples)

  stop(
    "\nDuplicate SampleID values detected."
  )

}


# ------------------------------------------------------------
# 8. VALIDATE EXPERIMENTAL CONDITIONS
# ------------------------------------------------------------

observed_conditions <- sort(
  unique(
    metadata$Condition
  )
)

required_conditions <- c(
  reference_condition,
  experimental_condition
)

unexpected_conditions <- setdiff(
  observed_conditions,
  required_conditions
)

if (length(unexpected_conditions) > 0) {

  stop(
    "\nUnexpected experimental conditions detected:\n",
    paste(
      unexpected_conditions,
      collapse = ", "
    ),
    "\n\nExpected only: ",
    paste(
      required_conditions,
      collapse = ", "
    )
  )

}

if (!all(
  required_conditions %in%
    observed_conditions
)) {

  stop(
    "\nOne or more required experimental conditions are absent."
  )

}


# ------------------------------------------------------------
# 9. STANDARDIZE METADATA FACTORS
# ------------------------------------------------------------

metadata <- metadata %>%

  mutate(

    CellLine = factor(
      CellLine
    ),

    Condition = factor(
      Condition,
      levels = c(
        reference_condition,
        experimental_condition
      )
    )

  )


# ------------------------------------------------------------
# 10. CHECK BIOLOGICAL REPLICATE STRUCTURE
# ------------------------------------------------------------

replicate_check <- metadata %>%

  group_by(
    CellLine,
    Condition
  ) %>%

  summarise(

    Biological_Replicates =
      n_distinct(
        BiologicalReplicate
      ),

    Samples =
      n(),

    .groups = "drop"

  )

print(
  replicate_check
)

if (
  any(
    replicate_check$Biological_Replicates != 3
  )
) {

  stop(
    "\nExpected exactly 3 biological replicates per ",
    "cell-line × condition group."
  )

}


# ------------------------------------------------------------
# 11. STOP SAFELY IF COUNT MATRIX IS NOT CONFIGURED
# ------------------------------------------------------------

if (is.null(count_matrix_file)) {

  message("\n============================================================")
  message("COUNT MATRIX NOT CONFIGURED")
  message("============================================================")

  message(
    "\nThe DESeq2 workflow is fully defined, but the count ",
    "matrix path has not yet been specified."
  )

  message(
    "\nSet the following object to the actual biological-level ",
    "count matrix:"
  )

  message(
    '\ncount_matrix_file <- "path/to/count_matrix.csv"'
  )

  message(
    "\nNo differential-expression results have been generated."
  )

  quit(
    save = "no",
    status = 0
  )
}


# ------------------------------------------------------------
# 12. VERIFY COUNT MATRIX FILE
# ------------------------------------------------------------

if (!file.exists(count_matrix_file)) {

  stop(
    "\nCount matrix file does not exist:\n",
    count_matrix_file
  )

}


# ------------------------------------------------------------
# 13. LOAD COUNT MATRIX
# ------------------------------------------------------------

message("\n============================================================")
message("Loading count matrix")
message("============================================================")

counts_df <- read_csv(
  count_matrix_file,
  show_col_types = FALSE
)

if (
  nrow(counts_df) == 0 ||
  ncol(counts_df) < 2
) {

  stop(
    "\nCount matrix appears to be empty or malformed."
  )

}


# ------------------------------------------------------------
# 14. IDENTIFY GENE-ID COLUMN
# ------------------------------------------------------------

gene_id_candidates <- c(
  "gene",
  "Gene",
  "gene_id",
  "Gene_ID",
  "GeneID",
  "ENSEMBL",
  "Ensembl",
  "GeneSymbol"
)

detected_gene_columns <- intersect(
  gene_id_candidates,
  colnames(counts_df)
)

if (length(detected_gene_columns) > 0) {

  gene_id_column <- detected_gene_columns[1]

} else {

  gene_id_column <- colnames(
    counts_df
  )[1]

  message(
    "\nNo standard gene-ID column detected."
  )

  message(
    "Using the first column as gene identifier: ",
    gene_id_column
  )

}


# ------------------------------------------------------------
# 15. EXTRACT GENE IDENTIFIERS
# ------------------------------------------------------------

gene_ids <- counts_df[
  [gene_id_column]
]

if (any(is.na(gene_ids))) {

  stop(
    "\nNA gene identifiers detected."
  )

}

gene_ids <- as.character(
  gene_ids
)

if (any(
  trimws(gene_ids) == ""
)) {

  stop(
    "\nEmpty gene identifiers detected."
  )

}


# ------------------------------------------------------------
# 16. HANDLE DUPLICATE GENE IDENTIFIERS
# ------------------------------------------------------------
#
# DESeq2 requires unique row identifiers.
#
# Duplicated identifiers are not silently discarded.
# They are aggregated by summing counts across rows.
#
# This is appropriate when duplicate rows represent the same
# gene identifier after upstream annotation.
# ------------------------------------------------------------

duplicate_gene_count <- sum(
  duplicated(gene_ids)
)

if (duplicate_gene_count > 0) {

  message(
    "\nDuplicate gene identifiers detected: ",
    duplicate_gene_count
  )

  message(
    "Duplicate rows will be aggregated by summing counts."
  )

}


# ------------------------------------------------------------
# 17. BUILD NUMERIC COUNT MATRIX
# ------------------------------------------------------------

count_data <- counts_df %>%

  select(
    -all_of(gene_id_column)
  ) %>%

  as.data.frame()


# Check numeric columns
non_numeric_columns <- names(
  count_data
)[
  !vapply(
    count_data,
    is.numeric,
    logical(1)
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

rownames(count_data) <- gene_ids


# ------------------------------------------------------------
# 18. AGGREGATE DUPLICATED GENES
# ------------------------------------------------------------

if (anyDuplicated(
  rownames(count_data)
) > 0) {

  count_data <- rowsum(
    count_data,
    group = rownames(
      count_data
    ),
    reorder = FALSE
  )

}


# ------------------------------------------------------------
# 19. VALIDATE COUNT VALUES
# ------------------------------------------------------------

count_matrix_numeric <- as.matrix(
  count_data
)

if (any(
  is.na(count_matrix_numeric)
)) {

  stop(
    "\nNA values detected in count matrix."
  )

}

if (any(
  !is.finite(count_matrix_numeric)
)) {

  stop(
    "\nNon-finite values detected in count matrix."
  )

}

if (any(
  count_matrix_numeric < 0
)) {

  stop(
    "\nNegative count values detected."
  )

}


# ------------------------------------------------------------
# 20. INTEGER CHECK
# ------------------------------------------------------------

non_integer_values <- sum(
  count_matrix_numeric !=
    round(count_matrix_numeric)
)

if (non_integer_values > 0) {

  warning(
    "\nThe count matrix contains ",
    non_integer_values,
    " non-integer entries."
  )

  message(
    "Values will be rounded before DESeq2 analysis."
  )

}

count_data <- round(
  count_data
)


# ------------------------------------------------------------
# 21. SAMPLE MATCHING
# ------------------------------------------------------------

metadata_samples <- metadata$SampleID

count_samples <- colnames(
  count_data
)

missing_samples <- setdiff(
  metadata_samples,
  count_samples
)

extra_samples <- setdiff(
  count_samples,
  metadata_samples
)

if (length(missing_samples) > 0) {

  stop(
    "\nSamples present in metadata but missing from count matrix:\n",
    paste(
      missing_samples,
      collapse = ", "
    )
  )

}

if (length(extra_samples) > 0) {

  warning(
    "\nAdditional count-matrix columns not present in metadata:\n",
    paste(
      extra_samples,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 22. RESTRICT TO ANALYSIS SAMPLES
# ------------------------------------------------------------

count_data <- count_data[
  ,
  metadata_samples,
  drop = FALSE
]


# ------------------------------------------------------------
# 23. FINAL SAMPLE-ORDER VALIDATION
# ------------------------------------------------------------

if (!identical(
  colnames(count_data),
  metadata_samples
)) {

  stop(
    "\nCount matrix and metadata sample order do not match."
  )

}

message(
  "✓ Count matrix successfully matched to metadata."
)


# ------------------------------------------------------------
# 24. LOW-COUNT FILTERING
# ------------------------------------------------------------
#
# Keep genes with >=10 counts in at least 3 biological samples.
#
# This removes genes with extremely limited information while
# retaining genes expressed in at least one biological group.
# ------------------------------------------------------------

minimum_count <- 10
minimum_samples <- 3

keep <- rowSums(
  count_data >= minimum_count
) >= minimum_samples

filtered_counts <- count_data[
  keep,
  ,
  drop = FALSE
]

filter_summary <- tibble(

  Genes_Before =
    nrow(count_data),

  Genes_After =
    nrow(filtered_counts),

  Genes_Removed =
    nrow(count_data) -
    nrow(filtered_counts),

  Minimum_Count =
    minimum_count,

  Minimum_Samples =
    minimum_samples

)

write_csv(
  filter_summary,
  file.path(
    results_dir,
    "DE_filtering_summary.csv"
  )
)


# ------------------------------------------------------------
# 25. CREATE DESEQ2 DATASET
# ------------------------------------------------------------

col_data <- metadata %>%

  select(
    SampleID,
    CellLine,
    Condition,
    BiologicalReplicate
  ) %>%

  column_to_rownames(
    "SampleID"
  )

dds <- DESeqDataSetFromMatrix(

  countData =
    filtered_counts,

  colData =
    col_data,

  design =
    ~ CellLine + Condition

)


# ------------------------------------------------------------
# 26. SET REFERENCE LEVELS
# ------------------------------------------------------------

dds$Condition <- relevel(
  dds$Condition,
  ref = reference_condition
)

dds$CellLine <- droplevels(
  dds$CellLine
)


# ------------------------------------------------------------
# 27. RUN DESEQ2
# ------------------------------------------------------------

message("\n============================================================")
message("Running DESeq2")
message("============================================================")

dds <- DESeq(
  dds,
  quiet = FALSE
)


# ------------------------------------------------------------
# 28. GLOBAL 3D VS 2D ANALYSIS
# ------------------------------------------------------------
#
# Model:
#
#   ~ CellLine + Condition
#
# The Condition coefficient therefore estimates the 3D versus
# 2D effect while accounting for systematic differences between
# the four cell lines.
# ------------------------------------------------------------

message("\n============================================================")
message("Global differential expression: 3D vs 2D")
message("============================================================")

global_results <- results(
  dds,
  contrast = c(
    "Condition",
    experimental_condition,
    reference_condition
  ),
  alpha = alpha
)

global_results_df <- as.data.frame(
  global_results
) %>%
  rownames_to_column(
    "Gene"
  ) %>%
  arrange(
    padj
  )


# ------------------------------------------------------------
# 29. ADD SIGNIFICANCE CATEGORIES
# ------------------------------------------------------------

global_results_df <- global_results_df %>%

  mutate(

    Significant_FDR =
      !is.na(padj) &
      padj < alpha,

    Significant_Effect =
      !is.na(log2FoldChange) &
      abs(log2FoldChange) >=
        lfc_threshold,

    Significant =
      Significant_FDR &
      Significant_Effect,

    Direction = case_when(

      Significant &
        log2FoldChange > 0 ~
        "Up_in_3D",

      Significant &
        log2FoldChange < 0 ~
        "Down_in_3D",

      TRUE ~
        "Not_significant"

    )

  )


# ------------------------------------------------------------
# 30. SAVE GLOBAL RESULTS
# ------------------------------------------------------------

write_csv(
  global_results_df,
  file.path(
    results_dir,
    "global_3D_vs_2D.csv"
  )
)


# ------------------------------------------------------------
# 31. GLOBAL SUMMARY
# ------------------------------------------------------------

global_summary <- tibble(

  Genes_Tested =
    sum(
      !is.na(
        global_results_df$padj
      )
    ),

  FDR_Significant =
    sum(
      global_results_df$Significant_FDR,
      na.rm = TRUE
    ),

  Effect_Size_Significant =
    sum(
      global_results_df$Significant_Effect,
      na.rm = TRUE
    ),

  Significant_Both =
    sum(
      global_results_df$Significant,
      na.rm = TRUE
    ),

  Up_in_3D =
    sum(
      global_results_df$Direction ==
        "Up_in_3D"
    ),

  Down_in_3D =
    sum(
      global_results_df$Direction ==
        "Down_in_3D"
    )

)

print(
  global_summary
)

write_csv(
  global_summary,
  file.path(
    results_dir,
    "global_DE_summary.csv"
  )
)


# ------------------------------------------------------------
# 32. TOP GLOBAL GENES
# ------------------------------------------------------------

top_global_genes <- global_results_df %>%

  filter(
    !is.na(padj)
  ) %>%

  arrange(
    padj
  ) %>%

  slice_head(
    n = 50
  )

write_csv(
  top_global_genes,
  file.path(
    results_dir,
    "top_50_global_DE_genes.csv"
  )
)


# ------------------------------------------------------------
# 33. WITHIN-CELL-LINE ANALYSIS FUNCTION
# ------------------------------------------------------------
#
# Each cell line is analyzed independently:
#
#   3D vs 2D
#
# This avoids assuming that the magnitude of the 3D response
# is identical across cell lines.
# ------------------------------------------------------------

run_cellline_DE <- function(
    count_matrix,
    sample_metadata,
    cell_line_name,
    reference = "2D",
    treatment = "3D",
    alpha = 0.05,
    lfc_threshold = 1
) {

  message(
    "\n------------------------------------------------------------"
  )

  message(
    "Cell line: ",
    cell_line_name
  )

  message(
    "------------------------------------------------------------"
  )


  # Subset metadata
  meta_sub <- sample_metadata %>%

    filter(
      CellLine == cell_line_name
    ) %>%

    droplevels()


  # Verify both conditions exist
  if (
    !all(
      c(reference, treatment) %in%
        meta_sub$Condition
    )
  ) {

    stop(
      "Both 2D and 3D conditions are required for ",
      cell_line_name,
      "."
    )

  }


  # Verify exactly 3 replicates per condition
  replicate_counts <- meta_sub %>%

    count(
      Condition,
      name = "n"
    )

  if (
    any(
      replicate_counts$n != 3
    )
  ) {

    stop(
      "Expected 3 biological replicates per condition for ",
      cell_line_name,
      "."
    )

  }


  # Sample IDs
  sample_ids <- meta_sub$SampleID


  # Subset count matrix
  counts_sub <- count_matrix[
    ,
    sample_ids,
    drop = FALSE
  ]


  # Metadata row order
  meta_sub <- meta_sub %>%

    select(
      SampleID,
      Condition,
      BiologicalReplicate
    ) %>%

    column_to_rownames(
      "SampleID"
    )


  # Ensure identical ordering
  counts_sub <- counts_sub[
    ,
    rownames(meta_sub),
    drop = FALSE
  ]


  # Create independent DESeq2 object
  dds_sub <- DESeqDataSetFromMatrix(

    countData =
      counts_sub,

    colData =
      meta_sub,

    design =
      ~ Condition

  )


  # Set reference
  dds_sub$Condition <- relevel(
    dds_sub$Condition,
    ref = reference
  )


  # Run DESeq2
  dds_sub <- DESeq(
    dds_sub,
    quiet = TRUE
  )


  # Extract contrast
  res_sub <- results(

    dds_sub,

    contrast = c(
      "Condition",
      treatment,
      reference
    ),

    alpha = alpha

  )


  # Convert to table
  res_df <- as.data.frame(
    res_sub
  ) %>%

    rownames_to_column(
      "Gene"
    ) %>%

    mutate(
      CellLine =
        cell_line_name,

      Significant_FDR =
        !is.na(padj) &
        padj < alpha,

      Significant_Effect =
        !is.na(log2FoldChange) &
        abs(log2FoldChange) >=
          lfc_threshold,

      Significant =
        Significant_FDR &
        Significant_Effect,

      Direction = case_when(

        Significant &
          log2FoldChange > 0 ~
          "Up_in_3D",

        Significant &
          log2FoldChange < 0 ~
          "Down_in_3D",

        TRUE ~
          "Not_significant"

      )

    ) %>%

    arrange(
      padj
    )


  return(
    res_df
  )
}


# ------------------------------------------------------------
# 34. RUN ALL CELL-LINE CONTRASTS
# ------------------------------------------------------------

cell_lines <- levels(
  metadata$CellLine
)

within_cellline_results <- map_dfr(

  cell_lines,

  ~ run_cellline_DE(

    count_matrix =
      filtered_counts,

    sample_metadata =
      metadata,

    cell_line_name =
      .x,

    reference =
      reference_condition,

    treatment =
      experimental_condition,

    alpha =
      alpha,

    lfc_threshold =
      lfc_threshold

  )

)


# ------------------------------------------------------------
# 35. SAVE WITHIN-CELL-LINE RESULTS
# ------------------------------------------------------------

write_csv(

  within_cellline_results,

  file.path(
    results_dir,
    "within_cell_line_3D_vs_2D.csv"
  )

)


# ------------------------------------------------------------
# 36. WITHIN-CELL-LINE SUMMARY
# ------------------------------------------------------------

within_summary <- within_cellline_results %>%

  group_by(
    CellLine
  ) %>%

  summarise(

    Genes_Tested =
      sum(
        !is.na(padj)
      ),

    FDR_Significant =
      sum(
        Significant_FDR,
        na.rm = TRUE
      ),

    Significant_Both =
      sum(
        Significant,
        na.rm = TRUE
      ),

    Up_in_3D =
      sum(
        Direction ==
          "Up_in_3D"
      ),

    Down_in_3D =
      sum(
        Direction ==
          "Down_in_3D"
      ),

    .groups = "drop"

  )

print(
  within_summary
)

write_csv(
  within_summary,
  file.path(
    results_dir,
    "within_cell_line_DE_summary.csv"
  )
)


# ------------------------------------------------------------
# 37. TOP GENES PER CELL LINE
# ------------------------------------------------------------

top_within_cellline <- within_cellline_results %>%

  filter(
    !is.na(padj)
  ) %>%

  group_by(
    CellLine
  ) %>%

  arrange(
    padj,
    .by_group = TRUE
  ) %>%

  slice_head(
    n = 25
  ) %>%

  ungroup()

write_csv(
  top_within_cellline,
  file.path(
    results_dir,
    "top_25_DE_genes_per_cell_line.csv"
  )
)


# ------------------------------------------------------------
# 38. SIGNIFICANT GENE TABLE
# ------------------------------------------------------------

significant_global <- global_results_df %>%

  filter(
    Significant
  ) %>%

  arrange(
    padj
  )

write_csv(
  significant_global,
  file.path(
    results_dir,
    "significant_global_DE_genes.csv"
  )
)


# ------------------------------------------------------------
# 39. UPREGULATED / DOWNREGULATED TABLES
# ------------------------------------------------------------

upregulated_global <- significant_global %>%

  filter(
    Direction == "Up_in_3D"
  )

downregulated_global <- significant_global %>%

  filter(
    Direction == "Down_in_3D"
  )

write_csv(
  upregulated_global,
  file.path(
    results_dir,
    "global_upregulated_in_3D.csv"
  )
)

write_csv(
  downregulated_global,
  file.path(
    results_dir,
    "global_downregulated_in_3D.csv"
  )
)


# ------------------------------------------------------------
# 40. SAVE DESEQ2 OBJECT
# ------------------------------------------------------------

saveRDS(
  dds,
  file = file.path(
    results_dir,
    "DESeq2_global_object.rds"
  )
)


# ------------------------------------------------------------
# 41. SAVE FILTERED COUNTS
# ------------------------------------------------------------

filtered_counts_export <- as.data.frame(
  filtered_counts
) %>%

  rownames_to_column(
    "Gene"
  )

write_csv(
  filtered_counts_export,
  file.path(
    results_dir,
    "filtered_counts_for_DESeq2.csv"
  )
)


# ------------------------------------------------------------
# 42. ANALYSIS PARAMETERS
# ------------------------------------------------------------

analysis_parameters <- tibble(

  Parameter = c(
    "Reference condition",
    "Experimental condition",
    "FDR threshold",
    "Absolute log2FC threshold",
    "Minimum count",
    "Minimum number of samples",
    "Global model",
    "Within-cell-line model"
  ),

  Value = c(
    reference_condition,
    experimental_condition,
    alpha,
    lfc_threshold,
    minimum_count,
    minimum_samples,
    "~ CellLine + Condition",
    "~ Condition"
  )

)

write_csv(
  analysis_parameters,
  file.path(
    results_dir,
    "DE_analysis_parameters.csv"
  )
)


# ------------------------------------------------------------
# 43. SESSION INFORMATION
# ------------------------------------------------------------

capture.output(

  sessionInfo(),

  file = file.path(
    results_dir,
    "DESeq2_sessionInfo.txt"
  )

)


# ------------------------------------------------------------
# 44. FINAL SUMMARY
# ------------------------------------------------------------

message("\n============================================================")
message("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED")
message("============================================================")

message(
  "\nGlobal 3D vs 2D:"
)

message(
  "  Genes tested          : ",
  global_summary$Genes_Tested
)

message(
  "  FDR < 0.05            : ",
  global_summary$FDR_Significant
)

message(
  "  FDR + |log2FC| >= 1  : ",
  global_summary$Significant_Both
)

message(
  "  Upregulated in 3D     : ",
  global_summary$Up_in_3D
)

message(
  "  Downregulated in 3D   : ",
  global_summary$Down_in_3D
)

message(
  "\nWithin-cell-line analyses completed for:"
)

for (
  cell_line in cell_lines
) {

  message(
    "  ✓ ",
    cell_line
  )

}

message(
  "\nSignificance definition:"
)

message(
  "  FDR < ",
  alpha,
  " AND |log2FC| >= ",
  lfc_threshold
)

message(
  "\nAnalysis model:"
)

message(
  "  Global        : ~ CellLine + Condition"
)

message(
  "  Cell-line     : ~ Condition"
)

message(
  "\nTechnical replicates were not treated as independent ",
  "biological replicates."
)

message("============================================================\n")
