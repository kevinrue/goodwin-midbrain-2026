```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR=/ceph/project/goodwin
RUN_DIR=$PROJECT_DIR/albrecht/reference_genome
```

Run the workflow:

```bash
cd $RUN_DIR
snakemake --dry-run --directory $RUN_DIR
nohup snakemake --directory $RUN_DIR > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```
