# Methodology

## Study design

This repository presents an exploratory transcriptomic analysis of colorectal cancer cell-line models comparing 3D and 2D growth conditions.

Four colorectal cancer cell lines were analyzed:

* HCT116
* HT29
* LS174T
* LS513

Each cell line included three biological replicates per growth condition. Technical replicates were collapsed to the biological-sample level before downstream differential expression analysis.

## Analysis workflow

The analysis followed the following workflow:

1. Dataset and sample identification
2. Metadata curation
3. Expression data preparation
4. Quality control
5. Technical replicate collapsing
6. Gene-level filtering
7. Differential expression analysis
8. Global 3D versus 2D comparison
9. Cell-line-specific 3D versus 2D comparisons
10. Chemokine–receptor panel analysis
11. Visualization and interpretation

## Quality control

Sample-level quality was evaluated using library-size summaries, principal component analysis (PCA), and sample-correlation analysis.

These analyses were used to assess overall sample structure, identify potential outliers, and examine whether biological replicates showed coherent expression patterns.

## Differential expression

Differential expression analysis was performed at the biological-sample level using DESeq2.

The primary comparison was:

**3D versus 2D**

Analyses were performed both globally across the four cell lines and separately within each cell line.

Genes were considered significant using a false discovery rate (FDR) threshold of 0.05 together with an absolute log2 fold-change threshold of 1.

## Chemokine–receptor panel

A predefined 18-gene chemokine–receptor panel was examined to characterize expression patterns across cell lines and growth conditions.

The panel included the following axes:

* CXCL9 / CXCL10 / CXCL11 – CXCR3
* CCL3 / CCL4 / CCL5 – CCR5
* CCL2 / CCL7 / CCL8 – CCR2
* CXCL12 – CXCR4
* CXCL16 – CXCR6
* CX3CL1 – CX3CR1

Panel-level analyses were performed using the biological-sample expression data and the differential-expression results.

## Interpretation

The analysis is intended as an exploratory computational study of colorectal cancer cell-line models.

The results describe transcriptional differences associated with 3D versus 2D growth conditions and should not be interpreted as clinical biomarker validation or direct evidence of patient-level biological effects.
