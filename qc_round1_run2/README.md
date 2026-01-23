```bash
conda activate snakemake
```

```bash
snakemake --dry-run
nohup snakemake --sdm conda &
```

```bash
deploy_root_path="/ceph/project/goodwin/datashare/albrecht/$(basename $(pwd))"
mkdir -p $deploy_root_path
rsync -avz reports/* $deploy_root_path
```
