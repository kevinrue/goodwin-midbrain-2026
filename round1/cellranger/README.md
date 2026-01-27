```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR="/ceph/project/goodwin"
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
WORKFLOW_DIR="$PROJECT_DIR/albrecht/workflow-cellranger_multi"
RUN_DIR="$PROJECT_DIR/albrecht/round1/cellranger"
export CELLRANGER_TARBALL="$PROJECT_DIR/albrecht/resources/cellranger-9.0.1.tar.gz"
```

Configure the workflow:

```bash
cd $RUN_DIR
mkdir data
# add symlinks to FASTQ files
cp -r $WORKFLOW_DIR/config $RUN_DIR
# edit as needed
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
  --set-resources cellranger_multi_run:runtime=720 \
  > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```
