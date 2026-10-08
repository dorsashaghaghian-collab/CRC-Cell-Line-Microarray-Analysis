# # CRC Cell-Line Microarray Analysis

### Exploratory transcriptomic analysis of colorectal cancer cell-line models using publicly available GEO datasets

---

## Overview

Colorectal cancer (CRC) is a biologically heterogeneous disease characterized by substantial variation in gene expression and cellular signaling programs.

This project focused on the exploratory analysis of publicly available transcriptomic datasets generated from colorectal cancer cell-line models. The primary goal was to evaluate the suitability of cell-line transcriptomic data for downstream investigation of gene-expression patterns and chemokine–receptor signaling in the context of CRC.

The project involved systematic dataset identification and screening, metadata inspection, sample classification, expression-data retrieval, preprocessing, and downstream transcriptomic analysis.

Rather than treating all publicly available datasets as equivalent, datasets were evaluated according to their experimental design, biological model, sample type, and relevance to the research question.

---

## Research Question

The initial objective was to determine whether publicly available colorectal cancer cell-line transcriptomic datasets could provide a suitable experimental framework for investigating gene-expression patterns associated with chemokine signaling.

Particular attention was given to chemokine–receptor axes of biological interest in colorectal cancer:

| Chemokine | Receptor |
|---|---|
| CXCL9 | CXCR3 |
| CXCL10 | CXCR3 |
| CXCL11 | CXCR3 |
| CCL3 | CCR5 |
| CCL4 | CCR5 |
| CCL5 | CCR5 |
| CCL2 | CCR2 |
| CCL7 | CCR2 |
| CCL8 | CCR2 |
| CXCL12 | CXCR4 |
| CXCL16 | CXCR6 |
| CX3CL1 | CX3CR1 |

These genes were subsequently considered as targeted biological features within the broader transcriptomic analysis.

---

## Project Objectives

The project was designed around several analytical objectives:

1. Identify publicly available GEO datasets containing relevant colorectal cancer cell-line transcriptomic data.
2. Screen candidate datasets according to sample type and experimental design.
3. Distinguish suitable cell-line datasets from datasets representing primary tissues, unrelated biological systems, or incompatible experimental designs.
4. Retrieve and organize relevant expression data and sample metadata.
5. Prepare the datasets for downstream transcriptomic analysis.
6. Investigate gene-expression patterns relevant to chemokine–receptor signaling.
7. Evaluate the suitability and limitations of cell-line models for the intended research question.

---

## Dataset Discovery and Screening

Candidate datasets were identified through the NCBI Gene Expression Omnibus (GEO).

Dataset selection was not based solely on the presence of colorectal cancer-related keywords. Candidate studies were examined with respect to:

- organism
- biological model
- cell-line identity
- experimental condition
- sample type
- sequencing or microarray platform
- availability of expression data
- sample metadata
- relevance to the research question

A total of multiple GEO studies were initially examined during the screening process, followed by progressive filtering according to experimental relevance and data availability.

### Dataset screening strategy

```text
GEO study identification
        ↓
Study-level metadata inspection
        ↓
Sample-level metadata inspection
        ↓
Cell-line eligibility assessment
        ↓
Expression-data availability
        ↓
Experimental-design assessment
        ↓
Candidate dataset selection
        ↓
Downstream transcriptomic analysis
