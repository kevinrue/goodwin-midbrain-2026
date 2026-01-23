stopifnot(require(BiocParallel))
stopifnot(require(dplyr))
stopifnot(require(DropletUtils))
stopifnot(require(readr))
stopifnot(require(stringr))
stopifnot(require(yaml))

if (interactive()) {
  workdir <- "/ceph/project/goodwin/albrecht/cellranger_round1_run2"
  raw_h5 <- file.path(workdir, "results/cellranger/WPPm024hrs_rep1/outs/multi/count/raw_feature_bc_matrix.h5")
  umi_stats_tsv <- "results/umi_stats/WPPm024hrs_rep1.tsv.gz"
  yaml <- "test.yaml"
  ncpus <- 1L
} else {
  raw_h5 <- snakemake@input[["raw"]]
  umi_stats_tsv <- snakemake@input[["umi_stats"]]
  yaml <- snakemake@output[["yaml"]]
  ncpus <- snakemake@threads
}

message("ncpus: ", ncpus)
if (identical(ncpus, 1L)) {
  BPPARAM <- SerialParam()
} else {
  BPPARAM <- MulticoreParam(
    workers = ncpus,
    jobname = "barcodeRanks",
    RNGseed = 1L
  )
}

sample_id <- basename(str_remove(raw_h5, "outs/multi/count/raw_feature_bc_matrix.h5"))

names(raw_h5) <- sample_id

sce <- read10xCounts(samples = raw_h5)

umi_stats <- read_tsv(umi_stats_tsv, show_col_types = FALSE)

keep_barcodes <- umi_stats |>
  filter(umi_sum > 0) |>
  pull(Barcode)

sce <- sce[, which(sce$Barcode %in% keep_barcodes)]

br_out <- barcodeRanks(
  m = sce,
  lower = 100,
  BPPARAM =
)


write_yaml(x = metadata(br_out), file = yaml)
