import os
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
import gzip
import shutil

# Function to generate standard compressed 10X format data for a single patient
def generate_standard_10x_all_gz(patient_name, output_dir, n_cells=500, n_genes=2000):
    np.random.seed(42)  # For reproducibility

    # Create fake gene names
    genes = [f"Gene_{i}" for i in range(n_genes)]
    marker_genes = ['EPCAM', 'CD3D', 'CD68', 'PECAM1']
    for i, gene in enumerate(marker_genes):
        genes[i] = gene

    # Create fake barcodes (cell IDs)
    barcodes = [f"{patient_name}_Cell_{i}" for i in range(n_cells)]

    # Generate counts matrix (sparse)
    data = np.random.poisson(0.5, size=(n_genes, n_cells))

    # Inject meaningful expression into marker genes
    for idx, marker in enumerate(marker_genes):
        if marker == 'EPCAM':
            data[idx, :int(0.3*n_cells)] += np.random.poisson(5, int(0.3*n_cells))
        elif marker == 'CD3D':
            data[idx, int(0.3*n_cells):int(0.6*n_cells)] += np.random.poisson(5, int(0.3*n_cells))
        elif marker == 'CD68':
            data[idx, int(0.6*n_cells):int(0.8*n_cells)] += np.random.poisson(5, int(0.2*n_cells))
        elif marker == 'PECAM1':
            data[idx, int(0.8*n_cells):] += np.random.poisson(5, int(0.2*n_cells))

    # Create output directory
    filtered_path = os.path.join(output_dir, patient_name, 'filtered_feature_bc_matrix')
    os.makedirs(filtered_path, exist_ok=True)

    # Save and compress matrix.mtx -> matrix.mtx.gz
    tmp_mtx_path = os.path.join(filtered_path, 'matrix.mtx')
    matrix = sparse.csr_matrix(data)
    mmwrite(tmp_mtx_path, matrix)

    with open(tmp_mtx_path, 'rb') as f_in, gzip.open(tmp_mtx_path + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_mtx_path)

    # Save and compress features.tsv -> features.tsv.gz
    tmp_feat_path = os.path.join(filtered_path, 'features.tsv')
    pd.DataFrame({
        'gene_id': genes,
        'gene_name': genes,
        'feature_type': ['Gene Expression'] * n_genes
    }).to_csv(tmp_feat_path, sep='\t', header=False, index=False)

    with open(tmp_feat_path, 'rb') as f_in, gzip.open(tmp_feat_path + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_feat_path)

    # Save and compress barcodes.tsv -> barcodes.tsv.gz
    tmp_barcode_path = os.path.join(filtered_path, 'barcodes.tsv')
    pd.DataFrame(barcodes).to_csv(tmp_barcode_path, sep='\t', header=False, index=False)

    with open(tmp_barcode_path, 'rb') as f_in, gzip.open(tmp_barcode_path + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_barcode_path)

    print(f"Standard 10X gzipped mock data generated for {patient_name} at {filtered_path}")

# Generate test datasets for multiple patients
for patient in ['Patient1', 'Patient2', 'Patient3']:
    generate_standard_10x_all_gz(patient_name=patient, output_dir='mock_scRNAseq_data', n_cells=500, n_genes=2000)


