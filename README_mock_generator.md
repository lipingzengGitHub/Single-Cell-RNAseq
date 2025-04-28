# scRNAseq_10x_mock_generator.py

This script generates **mock 10X Genomics single-cell RNA-seq data** for multiple samples,  
formatted exactly like real `filtered_feature_bc_matrix/` outputs (i.e., matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz).

It is designed for testing and benchmarking scRNA-seq analysis pipelines without using real datasets.

---

## Features

- Simulates multiple patient/sample datasets (e.g., Patient1, Patient2, Patient3).
- Generates 10X-compatible output:
  - `matrix.mtx.gz` (expression matrix in sparse Matrix Market format)
  - `features.tsv.gz` (gene identifiers and names)
  - `barcodes.tsv.gz` (cell barcodes)
- Includes basic biological realism by boosting expression for marker genes (e.g., EPCAM, CD3D, CD68, PECAM1) to simulate tumor cells, T cells, macrophages, and endothelial cells.
- Outputs ready for direct loading with `scanpy.read_10x_mtx()` or Seurat.

---

## Requirements

- Python 3.8+
- Libraries:
  ```bash
  pip install numpy pandas scipy

## Usage

python scRNAseq_10x_mock_generator.py

The output folder mock_scRNAseq_data/ will be created, containing: 

mock_scRNAseq_data/
├── Patient1/
│   └── filtered_feature_bc_matrix/
│       ├── matrix.mtx.gz
│       ├── features.tsv.gz
│       └── barcodes.tsv.gz
├── Patient2/
└── Patient3/




