#TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/h5ad/lung_cancer_all_oncoKB_19/"
#TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/h5ad/lung_cancer_oncoKB_with_drugs_11/"
#TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/h5ad/lung_cancer_all_oncoKB_mut_and_CGC_genes_30/"
TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/lung_cancer_all_oncoKB_19/"

import os
from anndata import read_h5ad
import pandas as pd

# for rerun lung cancer
TCGA_geneSubset_files = ["LungCancer_RNA_overlapped_genes_geneSubset.h5ad"]

for TCGA_file in TCGA_geneSubset_files:
    in_file = TCGA_path + TCGA_file
    out_file = TCGA_path + TCGA_file.replace(".h5ad", ".txt")
    adata = read_h5ad(in_file)
    df = adata.to_df().T
    df.to_csv(out_file, sep = "\t")
