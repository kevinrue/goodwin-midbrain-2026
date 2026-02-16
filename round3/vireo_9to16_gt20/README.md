```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR="/ceph/project/goodwin"
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
WORKFLOW_DIR="$PROJECT_DIR/albrecht/workflow-vireo"
RUN_DIR="$PROJECT_DIR/albrecht/round3/vireo_9to16_gt20"
DEPLOY_DIR="$PROJECT_DIR/datashare/albrecht/round3/vireo_9to16_gt20"
```

Configure the workflow:

```bash
cd $RUN_DIR
mkdir data
# add symlinks to FASTQ files
# edit config/* files as needed
```

Run the workflow:

```bash
cd $WORKFLOW_DIR
snakemake --dry-run \
  --sdm conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR
nohup snakemake \
  --sdm conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR \
  > $RUN_DIR/nohup.out 2>&1 &
tail -f $RUN_DIR/nohup.out
rm -rf $RUN_DIR/sps-*
```

Share the reports:

```bash
mkdir -p $DEPLOY_DIR
rsync -avz $RUN_DIR/reports/* $DEPLOY_DIR
```
