in_TCGA_path = "/home/yuching/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"
in_CCLE_path = "/home/yuching/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"

import anndata as ad
import os

all_files = os.listdir(in_TCGA_path)
RNA_files = []
for afile in all_files:
    if "_RNA.h5ad" in afile:
        RNA_files.append(afile)

for RNA_file in RNA_files:
    print(RNA_file)
    TCGA_data = ad.read_h5ad(in_TCGA_path + RNA_file)
    print(len(TCGA_data.var.index))
    CCLE_data = ad.read_h5ad(in_CCLE_path + RNA_file)
    print(len(CCLE_data.var.index))
    gene_list = list(set(CCLE_data.var.index).intersection(set(TCGA_data.var.index)))
    print(len(gene_list))
    TCGA_out_data = TCGA_data[:, gene_list]
    print(len(TCGA_out_data.var.index))
    CCLE_out_data = CCLE_data[:, gene_list]
    print(len(CCLE_out_data.var.index))
    out_RNA_file = RNA_file.replace("_RNA.h5ad", "_RNA_overlapped_genes.h5ad")
    TCGA_out_data.write_h5ad(in_TCGA_path + out_RNA_file)
    CCLE_out_data.write_h5ad(in_CCLE_path + out_RNA_file)
