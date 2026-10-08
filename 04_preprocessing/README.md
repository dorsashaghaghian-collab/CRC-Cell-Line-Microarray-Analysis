# 04 — Preprocessing

## Purpose

This stage documents the preprocessing procedures applied to transcriptomic expression data prior to downstream statistical and biological analysis.

The preprocessing workflow was designed to transform the acquired expression data into a consistent analytical format while preserving the biological information associated with each sample.

---

## Preprocessing Strategy

The preprocessing workflow included dataset-specific steps required to prepare the expression matrices for downstream analysis.

Conceptually:

```text
Raw / processed expression data
        ↓
Feature and sample inspection
        ↓
Identifier standardization
        ↓
Filtering / data cleaning
        ↓
Normalization or transformation
        ↓
Processed expression matrix
        ↓
Downstream analysis
