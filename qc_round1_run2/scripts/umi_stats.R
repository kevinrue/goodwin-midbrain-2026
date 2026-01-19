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
  gtf <- "/ceph/project/goodwin/albrecht/cellranger_mkref/results/cellranger_mkref/genes/genes.gtf.gz"
  tsv <- TODO
} else {
  raw_h5 <- snakemake@input[["raw"]]
  gtf <- snakemake@input[["gtf"]]
  tsv <- snakemake@output[["tsv"]]
}

sample_id <- basename(stringr::str_remove(raw_h5, "outs/multi/count/raw_feature_bc_matrix.h5"))

names(raw_h5) <- sample_id

gtf_data <- BiocIO::import(gtf, feature.type = "gene")

mt_gene_ids <- gtf_data$gene_id[which(seqnames(gtf_data) == "mitochondrion_genome")]

sce <- DropletUtils::read10xCounts(samples = raw_h5)

sce$umi_sum <- DelayedArray::colSums(counts(sce))

sce <- sce[, sce$umi_sum > 0]

sce$genes_detected <- DelayedArray::colSums(counts(sce) > 0)

sce$mt_pct <- ifelse(
  sce$umi_sum > 0,
  DelayedArray::colSums(counts(sce)[mt_gene_ids, ]) / sce$umi_sum,
  NA
)

write_tsv(
  as_tibble(colData(sce)),
  file = tsv
)
