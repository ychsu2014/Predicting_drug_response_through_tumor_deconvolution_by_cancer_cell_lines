import sys
#import scanpy as sc
from anndata import read_h5ad
import pandas as pd

in_path = sys.argv[1]
simu_file = sys.argv[2]
in_file = in_path + simu_file
out_file = in_path + simu_file[:-5] + ".txt"

# for raw testing data
#adata = sc.read_h5ad(in_file)
adata = read_h5ad(in_file)
df = adata.to_df().T
df.to_csv(out_file, sep = "\t")
