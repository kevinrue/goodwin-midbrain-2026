```bash
conda activate snakemake
```

```bash
snakedeploy deploy-workflow https://github.com/snakemake-workflows/cellranger-multi . --tag v2.0.0
```

```
export CELLRANGER_TARBALL="/ceph/project/goodwin/albrecht/resources/cellranger-9.0.1.tar.gz"
snakemake --dry-run
nohup snakemake --cores all --sdm conda --default-resources runtime=720 &
```
