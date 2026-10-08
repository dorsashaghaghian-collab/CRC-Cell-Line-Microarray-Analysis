# Analysis Report

## Overview

This project presents an exploratory transcriptomic analysis of colorectal cancer cell-line models grown under 3D and 2D conditions.

The analysis focused on four colorectal cancer cell lines:

* HCT116
* HT29
* LS174T
* LS513

The primary objective was to characterize transcriptional differences associated with 3D versus 2D growth and to examine a predefined chemokine–receptor gene panel within this experimental context.

---

## Experimental design

The final biological-level dataset contained:

| Feature                                         |   Design |
| ----------------------------------------------- | -------: |
| Cell lines                                      |        4 |
| Growth conditions                               |        2 |
| Biological replicates per cell line × condition |        3 |
| Final biological samples                        |       24 |
| Technical replicates per biological sample      |        8 |
| Chemokine–receptor panel                        | 18 genes |

Technical replicates were collapsed to the biological-sample level before differential expression analysis.

This resulted in a balanced design with three biological replicates for each cell-line × condition combination.

---

## Data processing and quality control

The initial expression matrix contained 58,243 genes.

After gene-level filtering, 19,515 genes remained for downstream analysis.

The minimum library size among the final biological samples was approximately 18.7 million reads.

Sample-level quality and structure were evaluated using:

* Library-size assessment
* Principal component analysis (PCA)
* Sample-correlation analysis

The PCA and correlation analyses were used to evaluate global sample structure and the consistency of biological replicates.

### PCA

The PCA provides a global view of transcriptional relationships among the biological samples.

Rather than interpreting individual genes from PCA, the analysis was used to assess whether samples showed coherent grouping according to experimental structure and whether any sample behaved as a potential outlier.

![Principal Component Analysis](../figures/PCA.png)

---

## Global differential expression analysis

A global comparison of 3D versus 2D growth conditions was performed using DESeq2.

A total of 19,474 genes were tested in the global differential-expression analysis.

Using:

* FDR < 0.05
* absolute log2 fold-change ≥ 1

the analysis identified:

**1,507 significant genes**

This indicates substantial transcriptional differences between the 3D and 2D growth conditions across the analyzed colorectal cancer cell-line models.

### Volcano plot

The volcano plot summarizes the global 3D-versus-2D differential-expression results by combining statistical significance and effect size.

![Global 3D versus 2D differential expression](../figures/volcano_3D_vs_2D.png)

---

## Cell-line-specific analysis

To determine whether transcriptional responses were shared across all cell lines or depended on cellular background, separate 3D-versus-2D comparisons were performed for each cell line.

| Cell line | Genes tested | Significant genes |
| --------- | -----------: | ----------------: |
| HCT116    |       17,702 |             1,100 |
| HT29      |       17,936 |             1,870 |
| LS174T    |       17,479 |               833 |
| LS513     |       16,693 |             1,033 |

The differences in the number of significant genes across cell lines indicate that the transcriptional response to 3D growth is not identical across the four models.

These results support the use of cell-line-specific analyses alongside the global comparison.

---

## Chemokine–receptor panel

A predefined 18-gene chemokine–receptor panel was examined across the biological samples.

The panel included:

* CXCL9 / CXCL10 / CXCL11 → CXCR3
* CCL3 / CCL4 / CCL5 → CCR5
* CCL2 / CCL7 / CCL8 → CCR2
* CXCL12 → CXCR4
* CXCL16 → CXCR6
* CX3CL1 → CX3CR1

The panel was evaluated independently against the biological-level expression matrix and then compared with the differential-expression results.

Importantly, absence of a gene from the final DESeq2-tested set was not interpreted as biological absence from the original dataset. Panel interpretation was therefore based on the available expression data and differential-expression results separately.

### Panel heatmap

The heatmap summarizes the expression landscape of the selected chemokine and receptor genes across the analyzed cell-line and growth-condition groups.

![Chemokine–receptor panel heatmap](../figures/chemokine_receptor_heatmap.png)

---

## Selected findings

Several genes in the chemokine–receptor panel showed significant changes in the global 3D-versus-2D comparison.

For example:

* **CXCL10** showed increased expression in 3D relative to 2D conditions (log2FC ≈ 1.40; FDR ≈ 0.046).
* **CCL5** showed a stronger 3D-associated increase (log2FC ≈ 1.93; FDR ≈ 3.7 × 10⁻⁹).

Cell-line-specific analyses also showed context-dependent responses.

For example, **CXCR4** showed a strong 3D-associated increase in HCT116, while CXCL10, CCL5, and CXCR4 displayed significant or condition-dependent changes in HT29 and LS513 in different comparisons.

These observations indicate that chemokine–receptor expression responses to 3D growth may depend on cellular background rather than representing a uniform response across all colorectal cancer models.

---

## Sample correlation

Sample-correlation analysis was used as an additional quality-control assessment of biological replicate consistency and overall sample relationships.

![Sample correlation](../figures/sample_correlation.png)

---

## Main conclusions

The analysis supports four main observations:

1. **3D growth is associated with substantial transcriptional remodeling** across colorectal cancer cell-line models.

2. **The magnitude of the transcriptional response differs between cell lines**, supporting cell-line-specific analysis rather than relying exclusively on a pooled comparison.

3. **Chemokine–receptor genes show condition- and cell-line-dependent expression patterns**, including significant changes involving CXCL10, CCL5, and CXCR4.

4. **The 18-gene panel provides a focused view of immune-signaling-related transcriptional changes**, but these findings remain exploratory and require validation in independent models and biological systems.

---

## Limitations

This study is an exploratory analysis of established colorectal cancer cell-line models.

Therefore:

* The findings cannot be directly interpreted as patient-level biomarkers.
* Cell-line-specific effects may reflect intrinsic prop
