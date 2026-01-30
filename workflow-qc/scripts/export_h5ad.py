import os
import anndata as ad
import pandas as pd
import scanpy as sc

ad.settings.allow_write_nullable_strings=True
# ad.settings.override(
#     allow_write_nullable_strings=True
# )

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

adata_filtered = adata[barcodes].copy()

# Add metadata
sample_id = snakemake.wildcards.sample
adata_filtered.obs['Sample'] = sample_id
adata_filtered.obs['Barcode'] = adata_filtered.obs.index.to_list()

# Write output
adata_filtered.write(output_h5ad)
