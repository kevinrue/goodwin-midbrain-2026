import pandas as pd

##### Helper functions #####

with open("/ceph/project/goodwin/albrecht/notebooks/gene_regions.bed") as fai:
    all_data = pd.read_table(
        fai,
        header=None
    )

# print(all_data)

# print(snakemake.wildcards['contig'])

filtered_data = all_data[all_data[0] == snakemake.wildcards['contig']]

# print(filtered_data)

filtered_data.to_csv(snakemake.output[0], sep='\t', index=False, header=False)
