# ============================================================
# CRC Cell-Line Transcriptomic Analysis
# Stage 01: GEO Dataset Screening
# ============================================================
#
# Purpose:
# Systematically screen candidate GEO studies to identify
# colorectal cancer cell-line transcriptomic datasets.
#
# Project context:
# This cell-line analysis was an exploratory branch of a broader
# colorectal cancer / chemokine-receptor investigation.
#
# Important:
# Cell-line datasets were evaluated as experimental models and
# were not treated as equivalent to primary human CRC tissue.
# ============================================================


# ------------------------------------------------------------
# 1. Required package
# ------------------------------------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("GEOquery")
}

library(GEOquery)


# ------------------------------------------------------------
# 2. Candidate GEO studies
# ------------------------------------------------------------
#
# These accessions represent studies examined during the
# dataset-screening process. Final eligibility was determined
# from study- and sample-level metadata.
# ------------------------------------------------------------

candidate_gse <- c(
  "GSE147456",
  "GSE159160",
  "GSE186791",
  "GSE193865",
  "GSE210045",
  "GSE226998",
  "GSE233326",
  "GSE252571",
  "GSE292302",
  "GSE309172",
  "GSE310994"
)


# ------------------------------------------------------------
# 3. Create output directory
# ------------------------------------------------------------

output_dir <- file.path("results")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 4. Function for GEO metadata retrieval
# ------------------------------------------------------------

get_geo_metadata <- function(gse_id) {

  message("Retrieving metadata for: ", gse_id)

  gse_object <- tryCatch(
    {
      getGEO(
        gse_id,
        GSEMatrix = TRUE,
        getGPL = FALSE
      )
    },
    error = function(e) {
      message("Could not retrieve ", gse_id)
      return(NULL)
    }
  )

  if (is.null(gse_object)) {
    return(NULL)
  }

  if (is.list(gse_object)) {
    eset <- gse_object[[1]]
  } else {
    eset <- gse_object
  }

  pheno <- pData(eset)

  list(
    gse = gse_id,
    eset = eset,
    pheno = pheno
  )
}


# ------------------------------------------------------------
# 5. Retrieve candidate datasets
# ------------------------------------------------------------

geo_results <- lapply(
  candidate_gse,
  get_geo_metadata
)

names(geo_results) <- candidate_gse


# ------------------------------------------------------------
# 6. Summarize study-level information
# ------------------------------------------------------------

screening_summary <- data.frame(
  GEO = character(),
  Samples = integer(),
  Platform = character(),
  Title = character(),
  stringsAsFactors = FALSE
)

for (gse_id in names(geo_results)) {

  result <- geo_results[[gse_id]]

  if (is.null(result)) {
    next
  }

  eset <- result$eset
  pheno <- result$pheno

  title_value <- ""

  if ("title" %in% colnames(pheno)) {
    title_value <- paste(
      unique(as.character(pheno$title)),
      collapse = " | "
    )
  }

  platform_value <- ""

  if ("platform_id" %in% colnames(pheno)) {
    platform_value <- paste(
      unique(as.character(pheno$platform_id)),
      collapse = "; "
    )
  }

  screening_summary <- rbind(
    screening_summary,
    data.frame(
      GEO = gse_id,
      Samples = nrow(pheno),
      Platform = platform_value,
      Title = title_value,
      stringsAsFactors = FALSE
    )
  )
}


# ------------------------------------------------------------
# 7. Inspect sample-level metadata
# ------------------------------------------------------------

metadata_columns <- lapply(
  geo_results,
  function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    colnames(x$pheno)
  }
)

metadata_columns


# ------------------------------------------------------------
# 8. Search metadata for cell-line related information
# ------------------------------------------------------------

cell_line_keywords <- paste(
  c(
    "cell line",
    "cell-line",
    "cellline",
    "colorectal",
    "colon",
    "CRC"
  ),
  collapse = "|"
)

cell_line_screening <- data.frame(
  GEO = character(),
  Sample = character(),
  MatchingMetadata = character(),
  stringsAsFactors = FALSE
)

for (gse_id in names(geo_results)) {

  result <- geo_results[[gse_id]]

  if (is.null(result)) {
    next
  }

  pheno <- result$pheno

  for (i in seq_len(nrow(pheno))) {

    sample_text <- paste(
      as.character(pheno[i, ]),
      collapse = " | "
    )

    if (
      grepl(
        cell_line_keywords,
        sample_text,
        ignore.case = TRUE
      )
    ) {

      cell_line_screening <- rbind(
        cell_line_screening,
        data.frame(
          GEO = gse_id,
          Sample = rownames(pheno)[i],
          MatchingMetadata = sample_text,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}


# ------------------------------------------------------------
# 9. Save screening outputs
# ------------------------------------------------------------

write.csv(
  screening_summary,
  file.path(
    output_dir,
    "GEO_study_screening_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_line_screening,
  file.path(
    output_dir,
    "cell_line_metadata_screening.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Display summary
# ------------------------------------------------------------

print(screening_summary)

cat("\n--------------------------------------------\n")
cat("Cell-line screening completed.\n")
cat("Studies evaluated:", nrow(screening_summary), "\n")
cat("Cell-line-related samples identified:",
    nrow(cell_line_screening), "\n")
cat("--------------------------------------------\n")


# ============================================================
# End of script
# ============================================================
