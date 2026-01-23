```bash
conda activate snakemake
```

```bash
WORKFLOW_DIR=/ceph/project/goodwin/albrecht/workflow-qc
RUN_DIR=/ceph/project/goodwin/albrecht/qc_round2
DEPLOY_DIR=/ceph/project/goodwin/datashare/albrecht/$(basename $RUN_DIR)
```

```bash
cd $WORKFLOW_DIR
snakemake --dry-run --directory $RUN_DIR
nohup snakemake --sdm conda --directory $RUN_DIR --conda-prefix ../conda > $RUN_DIR/nohup.out 2>&1 &
```

```bash
mkdir -p $DEPLOY_DIR
rsync -avz $RUN_DIR/reports/* $DEPLOY_DIR
```
