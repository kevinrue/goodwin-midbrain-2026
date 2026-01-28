```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR=/ceph/project/goodwin
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
WORKFLOW_DIR=$PROJECT_DIR/albrecht/workflow-qc
RUN_BASEDIR=round2/qc
RUN_DIR=$PROJECT_DIR/albrecht/$RUN_BASEDIR
DEPLOY_DIR=$PROJECT_DIR/datashare/albrecht/$RUN_BASEDIR
```

Run the workflow:

```bash
cd $WORKFLOW_DIR
snakemake --dry-run \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR
nohup snakemake \
  --sdm conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR \
  > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```

Share the reports:

```bash
mkdir -p $DEPLOY_DIR
rsync -avz $RUN_DIR/reports/* $DEPLOY_DIR
```
