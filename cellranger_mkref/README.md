```bash
conda activate snakemake
```

```bash
module load cellranger/9.0.0 # no snakemake workflow 
snakemake --dry-run
nohup snakemake --cores all --sdm conda apptainer &
```
