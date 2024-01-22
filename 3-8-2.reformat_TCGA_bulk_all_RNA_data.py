# from h5ad to txt
# the output txt file is for training/testing dataset preprocessing
in_path = "/home/yuching/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"

import os
from anndata import read_h5ad

all_files = os.listdir(in_path)
h5ad_files = []
for afile in all_files:
    if "_RNA_overlapped_genes.h5ad" in afile:
        h5ad_files.append(afile)

for afile in h5ad_files:
    in_file = in_path + afile
    adata = read_h5ad(in_file)
    df = adata.to_df().T
    out_file = in_file.replace(".h5ad", ".txt")
    df.to_csv(out_file, sep = "\t")
