import os
import pandas as pd
import scanpy as sc

# Get input/output from Snakemake
input_cellranger_h5 = snakemake.input.cellranger_h5
input_barcodes = snakemake.input.filtered_barcodes
output_h5ad = snakemake.output.h5ad

# Read cellranger H5
adata = sc.read_10x_h5(filename=input_cellranger_h5)

# Read and filter barcodes
barcodes = pd.read_csv(
    input_barcodes,
    header=None,
    names=['barcode']
)['barcode'].to_list()

adata = adata[barcodes].copy()

# Add metadata
sample_id = snakemake.wildcards.sample
adata.obs['Sample'] = sample_id
adata.obs['Barcode'] = adata.obs.index.to_list()

# Write output
adata.write(output_h5ad)
