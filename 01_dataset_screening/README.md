# 01 — Dataset Screening

## Purpose

This stage documents the identification and systematic screening of publicly available transcriptomic datasets relevant to colorectal cancer cell-line analysis.

The purpose of dataset screening was to identify studies with sufficient biological and technical compatibility for downstream transcriptomic analysis, rather than selecting datasets solely on the basis of disease-related keywords.

---

## Screening Criteria

Candidate GEO studies were evaluated according to the following criteria:

- Disease relevance: colorectal cancer or directly relevant experimental model
- Organism: primarily human datasets where appropriate
- Biological model: colorectal cancer cell lines
- Sample type and experimental context
- Availability of transcriptomic expression data
- RNA-seq or other compatible transcriptomic platforms
- Availability and quality of sample metadata
- Experimental group structure
- Suitability for downstream comparative analysis

---

## Screening Workflow

```text
GEO search
    ↓
Candidate study identification
    ↓
Study-level metadata inspection
    ↓
Sample-level metadata inspection
    ↓
Cell-line verification
    ↓
Experimental-design assessment
    ↓
Expression-data availability
    ↓
Eligibility assessment
    ↓
Candidate dataset selection
