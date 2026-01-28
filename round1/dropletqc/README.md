```bash
conda activate snakemake
```

Set up environment variables:

```bash
PROJECT_DIR="/ceph/project/goodwin"
PROJECT_CONDA_DIR="$PROJECT_DIR/albrecht/conda"
WORKFLOW_DIR="$PROJECT_DIR/albrecht/workflow-dropletqc"
RUN_DIR="$PROJECT_DIR/albrecht/round1/dropletqc"
export APPTAINER_CACHEDIR="/project/sims-lab/albrecht/mycontainers"
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
  --sdm apptainer conda \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR
nohup snakemake \
  --sdm apptainer conda \
  --apptainer-args "--bind $PROJECT_DIR" \
  --directory $RUN_DIR \
  --conda-prefix $PROJECT_CONDA_DIR \
  --set-resources dropletqc:mem=64G dropletqc:runtime=1h \
  > $RUN_DIR/nohup.out 2>&1 &
rm -rf $RUN_DIR/sps-*
```
