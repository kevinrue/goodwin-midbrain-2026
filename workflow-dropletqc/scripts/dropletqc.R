.libPaths(snakemake@input[["rlib"]])

stopifnot(require(DropletQC))
stopifnot(require(dplyr))
stopifnot(require(readr))
stopifnot(require(tibble))

input_bam <- snakemake@input[["bam"]]
input_barcodes <- snakemake@input[["barcodes"]]
params_id <- snakemake@params[["id"]]
output_tsv <- snakemake@output[["tsv"]]

message("input [bam]: ", input_bam)
message("input [barcodes]: ", input_barcodes)
message("params [id]: ", params_id)
message("output [tsv]: ", output_tsv)

nf_data <- nuclear_fraction_tags(
  bam = input_bam,
  barcodes = input_barcodes,
  tiles = 100,
  cores = snakemake@threads,
  verbose = FALSE
) |>
  rownames_to_column("Barcode") |>
  mutate(
    Sample = params_id,
    .before = 1
  )

write_tsv(
  nf_data,
  output_tsv
)
