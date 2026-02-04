```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR=/ceph/project/goodwin
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
RUN_DIR=$PROJECT_DIR/albrecht/reference_genome
```

Run the workflow:

```bash
cd $RUN_DIR
snakemake --dry-run \
  --sdm conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR
nohup snakemake \
  --sdm conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR \
  > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```
