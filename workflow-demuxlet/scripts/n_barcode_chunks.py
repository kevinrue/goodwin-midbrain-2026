import math
import pandas as pd

def n_barcode_chunks(file, size):
  print(f"requested chunks of size {size}")
  df = pd.read_csv(
    file,
    sep="\t",
    header=None,      # no header
    compression="gzip"
  )
  n_barcodes = len(df)
  print(f"{n_barcodes} barcodes in file")
  n_files = math.ceil(n_barcodes / size)
  print(f"will make {n_files} files")
  return n_files

n_barcode_chunks("/ceph/project/goodwin/albrecht/round2/qc/results/filter_barcodes/S1.tsv.gz", 1000)

