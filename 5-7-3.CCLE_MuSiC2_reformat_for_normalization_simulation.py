in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/h5ad/"

import os
import pandas as pd
import scanpy as sc
import numpy as np

all_files = os.listdir(in_path)
bulk_RNA_files = []
for afile in all_files:
    if "_pseudo_bulk.txt" in afile:
        bulk_RNA_files.append(afile)

def replaceGeneName(inGene):
    outGene = inGene.replace(".", "-")
    return(outGene)

##### the output MuSiC2 simulation data from R, "-" in gene names -> "."
def convertToh5ad(inFile, inPath = in_path, outPath = out_path):
    outFile = outPath + inFile.replace("_pseudo_bulk.txt", "_pseudo_bulk.h5ad")
    inFile = inPath + inFile
    inPhe = inFile.replace("_pseudo_bulk.txt", "_sampled_proportion.txt")
    df = pd.read_csv(inFile, sep = "\t")
    df = np.transpose(df)
    genes = df.columns
    new_genes = list(map(replaceGeneName, genes))
    column_name_replace = {}
    for i in range(len(new_genes)):
        old_gene = genes[i]
        new_gene = new_genes[i]
        column_name_replace[old_gene] = new_gene
    df.rename(columns = column_name_replace, inplace = True)
    data = sc.AnnData(df)
    phe_df = pd.read_csv(inPhe, sep = "\t")
    data.obs = phe_df
    data.write_h5ad(outFile)

for afile in bulk_RNA_files:
    convertToh5ad(afile)
