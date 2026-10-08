# 02 — Data Acquisition

## Purpose

This stage documents the acquisition and organization of transcriptomic datasets selected during the GEO screening process.

Following dataset-level and sample-level evaluation, eligible studies were identified for downstream analysis. The objective of this stage was to obtain the appropriate expression data and corresponding sample metadata while preserving the original study structure and experimental information.

---

## Data Source

All datasets analyzed in this project were obtained from publicly available studies deposited in the:

**NCBI Gene Expression Omnibus (GEO)**

Each study is identified by its unique GEO accession number (`GSEXXXXX`).

The original GEO accession is retained throughout the analysis to ensure traceability between the processed data and the source study.

---

## Acquisition Strategy

For each selected study, the following information was collected where available:

- GEO accession number
- study title
- organism
- experimental system
- cell-line identity
- sample identifiers
- experimental condition
- control/reference group
- sequencing or expression platform
- processed expression matrix
- sample metadata

The acquisition process was designed to preserve the relationship between expression measurements and their corresponding biological samples.

---

## Data Organization

The acquired datasets were organized according to their GEO accession and experimental structure.

Conceptually, the data flow was:

```text
GEO study
    │
    ├── Study metadata
    │
    ├── Sample metadata
    │
    └── Expression data
            │
            ▼
    Data organization
            │
            ▼
    Sample annotation
            │
            ▼
    Downstream preprocessing
