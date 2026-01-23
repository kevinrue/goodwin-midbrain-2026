```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR=/ceph/project/goodwin
WORKFLOW_DIR=$PROJECT_DIR/albrecht/workflow-qc
RUN_DIR=$PROJECT_DIR/albrecht/qc_round1_run2
DEPLOY_DIR=$PROJECT_DIR/datashare/albrecht/$(basename $RUN_DIR)
```

Run the workflow:

```bash
cd $WORKFLOW_DIR
snakemake --dry-run --directory $RUN_DIR --conda-prefix ../conda
nohup snakemake --sdm conda --directory $RUN_DIR --conda-prefix ../conda > $RUN_DIR/nohup.out 2>&1 &
```

Share the reports:

```bash
mkdir -p $DEPLOY_DIR
rsync -avz $RUN_DIR/reports/* $DEPLOY_DIR
```
