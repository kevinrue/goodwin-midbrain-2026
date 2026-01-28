# Set CRAN repository and install required packages
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Create and set library path
dir.create(snakemake@output[["path"]])
.libPaths(snakemake@output[["path"]])

# Install packages in library
install.packages(
  c(
    "devtools",
     "dplyr",
      "readr",
       "tibble"
  ),
  Ncpus = snakemake@threads
)
devtools::install_github(
  "powellgenomicslab/DropletQC",
  Ncpus = snakemake@threads
)
