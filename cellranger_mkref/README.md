```bash
conda activate snakemake
module load cellranger/9.0.0
```

Set up environment variables:

```bash
PROJECT_DIR=/ceph/project/goodwin
RUN_DIR=$PROJECT_DIR/albrecht/cellranger_mkref
```

Run the workflow:

```bash
cd $RUN_DIR
snakemake --dry-run
nohup snakemake \
  --directory $RUN_DIR \
  > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```
