# Single-Cell Transcriptomic Dissection of Human Invasive Ductal Carcinoma

## Overview
This project presents a complete single-cell RNA sequencing (scRNA-seq) analysis of 750 sorted cells from human Invasive Ductal Carcinoma (IDC) generated using 10x Genomics 3′ v3.1 chemistry and processed with Cell Ranger v6.0.0.

The workflow performs quality control, normalization using **SCTransform**, dimensionality reduction, clustering, automatic cell-type annotation, and tumor vs normal epithelial comparison to dissect the cellular heterogeneity of IDC.

---

## Dataset Information
- Tissue: Human Invasive Ductal Carcinoma (Breast Cancer)
- Cells: 750 sorted cells
- Platform:10x Genomics Chromium
- Chemistry: 3′ LT v3.1 (Universal 3′)
- Preprocessing: Cell Ranger 6.0.0

---

## Analysis Workflow
1. Data import from Cell Ranger output  
2. Quality control and filtering  
3. Normalization using SCTransform
4. PCA, graph-based clustering, and UMAP visualization  
5. Automatic cell-type annotation using SingleR 
6. Identification of epithelial cells  
7. Tumor vs normal epithelial state classification  
8. Differential gene expression analysis  
9. Visualization of marker genes and tumor-associated signatures  

---

## Tools & Packages
- Seurat
- SingleR
- celldex
- tidyverse
- patchwork

---

## Key Outputs
- UMAP plots of clustered and annotated cells  
- Cell-type annotations (immune, stromal, epithelial)  
- Tumor vs normal epithelial differential expression table  
- Marker gene heatmaps  
- Reproducible Seurat objects  

---

## Biological Insights
This analysis highlights:
- Intratumoral heterogeneity of IDC
- Distinct epithelial tumor and normal-like populations
- Immune and stromal components of the tumor microenvironment
- Tumor-associated transcriptional programs, including EMT-like states
