# Single-cell RNA-seq Analysis for Lung Cancer Patients

This project provides a complete scRNA-seq analysis pipeline for multiple lung cancer patient samples.  
It performs quality control, normalization, batch correction (Harmony), clustering, cell type annotation based on marker genes, differential expression analysis (DEGs) between tumor and non-tumor cells, and functional enrichment analysis (KEGG).

---

## Project Workflow

1. **Read scRNA-seq Data** (multiple patients)
2. **Quality Control** (filter low-quality cells and genes)
3. **Normalization & Log-transformation**
4. **Highly Variable Genes Selection**
5. **PCA for Dimensionality Reduction**
6. **Batch Correction using Harmony**
7. **Neighborhood Graph Construction and UMAP Visualization**
8. **Leiden Clustering**
9. **Cell Type Annotation** (e.g., Tumor cells, T cells, Macrophages, Endothelial cells)
10. **Differential Expression Analysis** (Tumor vs Other cell types)
11. **Functional Enrichment Analysis** (KEGG pathways)
12. **Visualization** (UMAP, Volcano Plot)

---

## Installation

You need Python 3.8+ and the following packages:

```bash
pip install scanpy harmonypy gseapy matplotlib seaborn

## Usage

Prepare your data:
Ensure you have 10X Genomics output directories (filtered_feature_bc_matrix) for each patient sample.

Update sample list:
Edit the sample names in the script:

samples = ['Patient1', 'Patient2', 'Patient3']

## How to Run
python Single-Cell-RNAseq_analysis.py

Note: Modify the data folder in the script

## Output

Tumor_vs_Other_DEGs_multiPatient.csv: Differentially expressed genes between tumor and non-tumor cells.

KEGG_enrichment/: KEGG pathway enrichment analysis results.

Volcano plots and UMAP visualizations.

