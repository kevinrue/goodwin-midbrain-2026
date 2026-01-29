```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR="/ceph/project/goodwin"
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
WORKFLOW_DIR="$PROJECT_DIR/albrecht/workflow-dna-seq-gatk-variant-calling"
RUN_DIR="$WORKFLOW_DIR"
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
rm -rf $RUN_DIR/sps-*
```
