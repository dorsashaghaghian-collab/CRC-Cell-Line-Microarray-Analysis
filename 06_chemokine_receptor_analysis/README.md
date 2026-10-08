# 06 — Chemokine–Receptor Analysis

## Purpose

This stage focuses on the targeted evaluation of chemokine and chemokine-receptor genes within the transcriptomic datasets analyzed in this project.

The analysis was designed to examine whether genes belonging to biologically relevant chemokine–receptor axes showed detectable and/or differential expression within colorectal cancer cell-line models.

---

## Biological Framework

The analysis focused on the following chemokine–receptor relationships:

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

These genes were selected as targeted features based on their relevance to chemokine-mediated cellular signaling.

---

## Analytical Strategy

The targeted analysis followed the general workflow:

```text
Processed expression matrix
        ↓
Gene identifier mapping
        ↓
Chemokine / receptor gene selection
        ↓
Expression extraction
        ↓
Group-wise comparison
        ↓
Differential-expression results
        ↓
Targeted visualization
        ↓
Biological interpretation
