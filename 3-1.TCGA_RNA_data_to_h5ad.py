in_RNA = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena_removeSLC35E2dup"
out_total = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/all_RNA.h5ad"
out_idx = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/all_RNA_idx_sID.txt"

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import anndata as ad

# duplicate entries for the "SLE35E2" gene, drop the one with lower  value
#>>> RNA.iloc[np.where(RNA["sample"] == "SLC35E2")]
#        sample  TCGA-OR-A5J1-01  TCGA-OR-A5J2-01  TCGA-OR-A5J3-01  TCGA-OR-A5J5-01  ...  TCGA-CG-4472-01  TCGA-CG-4474-01  TCGA-CG-4475-01  TCGA-CG-4476-01  TCGA-CG-4477-01
#16301  SLC35E2            11.69             9.70            10.28            10.19  ...            10.44            10.88            10.75            11.90            10.00
#16302  SLC35E2             5.18             4.48             3.97             3.01  ...             4.46             6.02             5.28             6.35             4.27

f = open(in_RNA)
fout = open(out_idx, "w")
lines = f.readlines()
header = lines[0].strip("\n").split("\t")
sID_list = header[1:]
sID_idx_dict = {}
count = -1
for sID in sID_list:
    count += 1
    sID_idx_dict[count] = sID

for i in sID_idx_dict.keys():
    sID = sID_idx_dict[i]
    fout.write(str(i) + "\t" + sID + "\n")
f.close()
fout.close()

lines = lines[1:]
gene_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    gene_list.append(gene)
gene_idx_dict = {}
count = -1
for gene in gene_list:
    count += 1
    gene_idx_dict[gene] = count

col_idx_list = []
row_idx_list = []
data_list = []
col_idx = -1
for line in lines:
    col_idx += 1
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    cols = cols[1:]
    row_idx = -1
    for i in range(len(cols)):
        row_idx += 1
        data = cols[i]
        col_idx_list.append(col_idx)
        row_idx_list.append(row_idx)
        if data != "NA":
            data_list.append(float(data))
        # let NA value be zero
        else:
            data_list.append(0)

row_idx_array = np.array(row_idx_list)
col_idx_array = np.array(col_idx_list)
adata = csr_matrix((data_list, (row_idx_array, col_idx_array)), shape = (max(row_idx_array)+1, max(col_idx_array)+1))
adata = ad.AnnData(adata, dtype = float)
adata.var_names = gene_list

# write all RNA data to file
adata.write_h5ad(out_total)