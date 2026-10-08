# CRC Cell-Line Transcriptomic Analysis

### Exploratory analysis of 3D versus 2D growth in colorectal cancer cell-line models

This repository contains an exploratory transcriptomic analysis of colorectal cancer cell-line models grown under **3D and 2D conditions**, with a focused analysis of chemokine–receptor expression patterns.

The analysis integrates sample-level quality control, differential expression analysis, cell-line-specific comparisons, and a predefined 18-gene chemokine–receptor panel.

---

## Research question

**How does 3D growth alter the transcriptional landscape of colorectal cancer cell-line models, and are chemokine–receptor expression patterns affected in a cell-line-dependent manner?**

---

## Study design

| Feature                  | Design                      |
| ------------------------ | --------------------------- |
| Dataset                  | GSE185055                   |
| Cell lines               | HCT116, HT29, LS174T, LS513 |
| Conditions               | 2D vs 3D                    |
| Biological replicates    | 3 per cell-line × condition |
| Final biological samples | 24                          |
| Technical replicates     | 8 per biological sample     |
| Focused panel            | 18 chemokine/receptor genes |
| Differential expression  | DESeq2                      |

Technical replicates were collapsed to the biological-sample level before downstream differential-expression analysis.

---

## Analysis workflow

```text
Dataset identification
        ↓
Metadata curation
        ↓
Expression data preparation
        ↓
Quality control
        ↓
Technical replicate collapsing
        ↓
Gene filtering
        ↓
Differential expression
        ↓
Global 3D vs 2D analysis
        ↓
Cell-line-specific analysis
        ↓
Chemokine–receptor panel
        ↓
Visualization and interpretation
```

---

# Key results

## Global transcriptional response

The global 3D-versus-2D comparison tested **19,474 genes**.

Using:

* FDR < 0.05
* |log2 fold-change| ≥ 1

the analysis identified **1,507 significant genes**.

This indicates substantial transcriptional remodeling associated with 3D growth across the analyzed colorectal cancer cell-line models.

### Global differential expression

![Global 3D versus 2D differential expression](figures/volcano_3D_vs_2D.png)

---

## Cell-line-specific response

The magnitude of the transcriptional response differed between cell lines:

| Cell line | Genes tested | Significant genes |
| --------- | -----------: | ----------------: |
| HCT116    |       17,702 |             1,100 |
| HT29      |       17,936 |             1,870 |
| LS174T    |       17,479 |               833 |
| LS513     |       16,693 |             1,033 |

This heterogeneity indicates that the transcriptional response to 3D growth is not uniform across the four colorectal cancer models.

---

# Chemokine–receptor analysis

A predefined **18-gene panel** was used to investigate selected chemokine–receptor axes:

| Receptor axis | Ligands               |
| ------------- | --------------------- |
| CXCR3         | CXCL9, CXCL10, CXCL11 |
| CCR5          | CCL3, CCL4, CCL5      |
| CCR2          | CCL2, CCL7, CCL8      |
| CXCR4         | CXCL12                |
| CXCR6         | CXCL16                |
| CX3CR1        | CX3CL1                |

### Chemokine–receptor expression landscape

![Chemokine–receptor panel heatmap](figures/chemokine_receptor_heatmap.png)

Several panel genes showed significant 3D-associated changes in the global analysis.

For example:

* **CXCL10:** log2FC ≈ 1.40, FDR ≈ 0.046
* **CCL5:** log2FC ≈ 1.93, FDR ≈ 3.7 × 10⁻⁹

Cell-line-specific analyses further indicated that some chemokine–receptor responses were dependent on cellular background.

For example, **CXCR4** showed a strong 3D-associated change in HCT116, while CXCL10, CCL5, and CXCR4 showed context-dependent changes across other cell lines.

---

# Quality control

## Principal component analysis

PCA was used to evaluate global sample structure and assess relationships among biological replicates.

![Principal Component Analysis](figures/PCA.png)

## Sample correlation

Sample-correlation analysis provided an additional assessment of replicate consistency and overall sample relationships.

![Sample correlation](figures/sample_correlation.png)

---

# Interpretation

The analysis supports several observations:

**1. 3D growth is associated with substantial transcriptional remodeling.**

More than 1,500 genes met the predefined significance and effect-size criteria in the global comparison.

**2. The response is cell-line dependent.**

The number of significant genes differed considerably among HCT116, HT29, LS174T, and LS513.

**3. Chemokine–receptor expression is context dependent.**

Selected chemokine and receptor genes showed 3D-associated changes, but these patterns were not identical across cell lines.

**4. The 18-gene panel provides a focused view of immune-signaling-related transcriptional changes.**

The panel is useful for exploratory characterization but does not constitute a clinical biomarker or mechanistic validation study.

---

# Repository structure

```text
CRC-Cell-Line-Microarray-Analysis/
│
├── README.md
│
├── figures/
│   ├── chemokine_receptor_heatmap.png
│   ├── PCA.png
│   ├── sample_correlation.png
│   └── volcano_3D_vs_2D.png
│
├── results/
│   ├── chemokine_receptor_panel.csv
│   ├── global_3D_vs_2D.csv
│   └── within_cell_line_3D_vs_2D.csv
│
├── data/
│   └── sample_metadata.csv
│
├── scripts/
│   ├── 01_metadata_and_design.R
│   ├── 02_quality_control.R
│   ├── 03_differential_expression.R
│   ├── 04_chemokine_receptor_analysis.R
│   └── 05_visualization.R
│
└── docs/
    ├── methodology.md
    └── analysis_report.md
```

---

# Reproducibility

The repository provides:

* Curated biological-level sample metadata
* Selected differential-expression results
* Chemokine–receptor panel results
* Quality-control and analysis figures
* R-based analysis workflow documentation
* Detailed methodology and analysis reports

Raw GEO data and large intermediate files are intentionally not stored in the repository.

The source dataset is publicly available through the **NCBI Gene Expression Omnibus (GEO), accession GSE185055**.

---

# Limitations

This project is an **exploratory computational analysis of established colorectal cancer cell-line models**.

Therefore:

* Results should not be interpreted as patient-level biomarker validation.
* Cell-line-specific effects may reflect intrinsic properties of individual models.
* The chemokine–receptor panel represents a predefined targeted analysis rather than a genome-wide mechanistic model.
* Statistical association does not establish causality.
* Independent datasets and experimental validation would be required for translational conclusions.

---

# Author

**Dorsa Shaghaghian**

B.Sc. Biotechnology, Yazd University

Research interests: computational biology, transcriptomics, cancer biology, functional genomics, and molecular biotechnology.

---

## Citation

If you use this repository or its analysis framework, please cite the repository and the underlying GEO dataset.

See `CITATION.cff` for repository citation information.
