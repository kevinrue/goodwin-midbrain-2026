stopifnot(require(dplyr))
stopifnot(require(readr))

input_stats <- snakemake@input[["umi_stats"]]
params_umi_min <- snakemake@params[["umi_min"]]
params_genes_min <- snakemake@params[["genes_min"]]
params_mt_max <- snakemake@params[["mt_max"]]
output_tsv <- snakemake@output[["tsv"]]

message("input [umi_stats]: ", input_stats)
message("params [umi_min]: ", params_umi_min)
message("params [genes_min]: ", params_genes_min)
message("params [mt_max]: ", params_mt_max)
message("output [tsv]: ", output_tsv)

raw_qc_data <- read_tsv(
  file = input_stats,
  show_col_types = FALSE
)
raw_qc_data |> 
  filter(
    cellranger_filtered == TRUE &
    umi_sum >= params_umi_min &
    genes_detected >= params_genes_min &
    mt_pct <= params_mt_max
  ) |>
  select(Barcode) |>
  write_tsv(
    file = output_tsv,
    col_names = FALSE
  )

