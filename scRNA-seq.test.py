#!/usr/bin/env python 3

#Install required tools
#pip install scanpy harmonypy gseapy matplotlib seaborn igraph leidenalg

import scanpy as sc
import pandas as pd
import numpy as np
import harmonypy as hm
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp

# Step 1. Read multiple patient samples
# Assume each patient sample is stored in a separate folder
samples = ['Patient1', 'Patient2', 'Patient3']
adatas = []

for sample in samples:
    ad = sc.read_10x_mtx(f'path/to/{sample}/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True) #Modify the path for your data
    ad.obs['batch'] = sample  # Add batch labels
    adatas.append(ad)

# Merge all samples into a single AnnData object
adata = adatas[0].concatenate(adatas[1:], batch_key="batch", batch_categories=samples)

# Step 2. Quality Control (QC)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Annotate mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Step 3. Filter good-quality cells
adata = adata[(adata.obs.n_genes_by_counts > 200) & (adata.obs.n_genes_by_counts < 6000) & (adata.obs.pct_counts_mt < 10), :]

# Step 4. Normalize total counts per cell and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Step 5. Identify highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# Step 6. Scale and perform initial PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# Step 7. Batch correction using Harmony
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_pca_harmony'] = ho.Z_corr.T

# Step 8. Build neighbors graph and run UMAP using Harmony-corrected PCs
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=10)
sc.tl.umap(adata)

# Step 9. Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

# Step 10. Plot markers to help annotate cell types
# Check whether marker genes existed or not
marker_genes = ['EPCAM', 'CD3D', 'CD68', 'PECAM1']
marker_genes_exist = [g for g in marker_genes if g in adata.var_names]
sc.pl.umap(adata, color=['leiden'] + marker_genes_exist)

#sc.pl.umap(adata, color=['leiden', 'EPCAM', 'CD3D', 'CD68', 'PECAM1'])

# Step 11. Enhanced cell type annotation based on multiple marker genes

# Define marker genes and thresholds
marker_genes_thresholds = {
    'Tumor': {'EPCAM': 0.5},
    'T_cell': {'CD3D': 0.5},
    'Macrophage': {'CD68': 0.5},
    'Endothelial': {'PECAM1': 0.5}
}

# Initialize
cell_types = []

# For each cell, determine type
for i in range(adata.n_obs):
    assigned_type = 'Other'
    for cell_type, markers in marker_genes_thresholds.items():
        passed = True
        for gene, threshold in markers.items():
            if gene in adata.var_names:
                if adata[i, gene].X.mean() <= threshold:
                    passed = False
                    break
            else:
                passed = False
                break
        if passed:
            assigned_type = cell_type
            break
    cell_types.append(assigned_type)

# Save result
adata.obs['cell_type'] = cell_types

# Optional: check summary
print(adata.obs['cell_type'].value_counts())


# Step 12. Differential expression analysis (Tumor vs Other cells)
sc.tl.rank_genes_groups(adata, 'cell_type', groups=['Tumor'], reference='Other', method='wilcoxon')
degs = sc.get.rank_genes_groups_df(adata, group='Tumor')
degs = degs[degs['pvals_adj'] < 0.05]
degs.to_csv('Tumor_vs_Other_DEGs_multiPatient.csv', index=False)

# Step 13. Functional enrichment analysis (KEGG)
gene_list = degs['names'].tolist()
enr = gp.enrichr(gene_list=gene_list, organism='Human', gene_sets=['KEGG_2021_Human'], outdir='KEGG_enrichment', cutoff=0.05)
gp.barplot(enr.res2d, title='KEGG Enrichment of Tumor DEGs')

# Step 14. Volcano plot visualization
degs['log10_pval'] = -np.log10(degs['pvals_adj'])
plt.figure(figsize=(8,6))
sns.scatterplot(data=degs, x='logfoldchanges', y='log10_pval', hue='pvals_adj', palette='coolwarm', edgecolor=None)
plt.title('Volcano Plot: Tumor vs Other Cells (multi-patient)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted p-value')
plt.grid(True)
plt.show()
























































