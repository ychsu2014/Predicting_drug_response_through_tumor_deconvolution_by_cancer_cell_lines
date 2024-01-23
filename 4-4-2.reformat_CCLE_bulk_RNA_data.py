in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"

import os
from anndata import read_h5ad
import pandas as pd

all_files = os.listdir(in_path)
in_geneSubset_files = []
for afile in all_files:
    if "_RNA_ovelapped_genes_geneSubset.h5ad" in afile:
        in_geneSubset_files.append(afile)

# for rerun skin cancer
#TCGA_geneSubset_files = ["SkinCancer_RNA_geneSubset.h5ad"]

for afile in in_geneSubset_files:
    in_file = in_path + afile
    out_file = in_file.replace(".h5ad", ".txt")
    adata = read_h5ad(in_file)
    df = adata.to_df().T
    df.to_csv(out_file, sep = "\t")
