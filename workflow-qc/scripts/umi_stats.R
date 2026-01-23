stopifnot(require(BiocIO))
stopifnot(require(DelayedArray))
stopifnot(require(DropletUtils))
stopifnot(require(rtracklayer))
stopifnot(require(readr))
stopifnot(require(stringr))
stopifnot(require(SingleCellExperiment))
stopifnot(require(SummarizedExperiment))
stopifnot(require(tibble))

if (interactive()) {
  workdir <- "/ceph/project/goodwin/albrecht/cellranger_round1_run2"
  raw_h5 <- file.path(workdir, "results/cellranger/WPPm024hrs_rep1/outs/multi/count/raw_feature_bc_matrix.h5")
  filtered_tsv <- file.path(workdir, "results/cellranger/WPPm024hrs_rep1/outs/per_sample_outs/WPPm024hrs_rep1/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz")
  gtf <- "/ceph/project/goodwin/albrecht/cellranger_mkref/results/cellranger_mkref/genes/genes.gtf.gz"
  tsv <- "test.tsv.gz"
} else {
  raw_h5 <- snakemake@input[["raw"]]
  filtered_tsv <- snakemake@input[["filtered"]]
  gtf <- snakemake@input[["gtf"]]
  tsv <- snakemake@output[["tsv"]]
}

sample_id <- basename(str_remove(raw_h5, "outs/multi/count/raw_feature_bc_matrix.h5"))

names(raw_h5) <- sample_id

gtf_data <- import(gtf, feature.type = "gene")

mt_gene_ids <- gtf_data$gene_id[which(seqnames(gtf_data) == "mitochondrion_genome")]

sce <- read10xCounts(samples = raw_h5)

sce$umi_sum <- colSums(counts(sce))

sce <- sce[, sce$umi_sum > 0]

sce$genes_detected <- colSums(counts(sce) > 0)

sce$mt_pct <- ifelse(
  sce$umi_sum > 0,
  colSums(counts(sce)[mt_gene_ids, ]) / sce$umi_sum,
  NA
)

cellranger_filtered_barcodes <- read_tsv(
  filtered_tsv,
  col_names = c("id"),
  show_col_types = FALSE
)

sce$cellranger_filtered <- sce$Barcode %in% cellranger_filtered_barcodes$id

write_tsv(
  colData(sce) |>
    as_tibble(),
  file = tsv
)
