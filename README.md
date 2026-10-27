# Single-Cell RNA-Seq Profiling of PBMC Microenvironments

End-to-end computational pipeline for processing, quality controlling, and annotating high-throughput single-cell RNA sequencing (scRNA-seq) data. Identified distinct immune cell populations from 3k Peripheral Blood Mononuclear Cells (PBMCs) using unsupervised machine learning and differential gene expression.

## Computational Pipeline
The analysis was conducted in R using the `Seurat` framework, adhering to standard genomic data science protocols:
1. **Quality Control:** Filtered dead/apoptotic cells (mitochondrial RNA > 5%) and multiplet droplets (feature count > 2500 or < 200).
2. **Normalization & Feature Selection:** Log-normalized data and extracted the top 2,000 highly variable genes (HVGs).
3. **Dimensionality Reduction:** 
   * **Linear:** Applied Principal Component Analysis (PCA) to evaluate feature loadings. Identified a strong myeloid vs. lymphoid variance signature on PC_1 (driven by LYZ/CST3).
   * **Non-Linear:** Projected the top 10 Principal Components into a 2D UMAP space for topological visualization.
4. **Graph-Based Clustering:** Constructed a K-Nearest Neighbor (KNN) graph and grouped cells using Louvain clustering.
5. **Biomarker Identification:** Extracted cluster-specific transcriptomic signatures using differential expression testing.

## Cell Type Annotation Results
Mathematical clusters were successfully mapped to known biological cell types using the following identified gene signatures:
* **Naïve CD4+ T Cells:** `CCR7`, `LEF1`
* **CD14+ Monocytes:** `S100A12`
* **B Cells:** `VPREB3`
* **Cytotoxic T/NK Cells:** `GZMK`, `GZMH`
* **Dendritic Cells:** `FCER1A`

## Repository Contents
* `scripts/01_seurat_pbmc_pipeline.R`: The complete R script containing the Seurat workflow.
* `plots/`: Exported visualizations including QC violin plots, PCA dimension loadings, and the final UMAP cell clusters.

## Technical Stack
* **Language:** R
* **Core Libraries:** `Seurat`, `dplyr`, `ggplot2`
