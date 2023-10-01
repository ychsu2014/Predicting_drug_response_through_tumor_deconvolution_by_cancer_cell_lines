# Scaden simulation method
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/"
# MuSiC2 simulation method
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/h5ad/"
in_gene_len = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/gene_length_table.txt"

import os
import anndata as ad
import pandas as pd
from bioinfokit.analys import norm
import numpy as np
from sklearn import preprocessing as pp
from scipy import sparse

all_files = os.listdir(in_path)
pseudo_bulk_files = []
for afile in all_files:
    if "_pseudo_bulk.h5ad" in afile:
        pseudo_bulk_files.append(afile)

gene_length_df = pd.read_csv(in_gene_len, sep = "\t", index_col = 0)

## for testing
#pseudo_bulk_files = ["Thyroid_Cancer_pseudo_bulk.h5ad"]

def sample_scaling(x):
    mms = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
    ## it scales features so transpose is needed
    x = mms.fit_transform(x.T.toarray()).T
    x = sparse.csr_matrix(x)
    return(x)
   
## to check if the calculated TPM values are correct:
#(temp_df.loc["A1BG"]["1"]/(temp_df.loc["A1BG"]["gene_length"]/1000))/((sum(temp_df["1"]/(temp_df["gene_length"]/1000)))/1000000)

for CCLE_pseudo_bulk_file in pseudo_bulk_files:
    print(CCLE_pseudo_bulk_file)
    out_h5ad = CCLE_pseudo_bulk_file.replace(".h5ad", "_norm_scaled.h5ad")
    data = ad.read_h5ad(in_path + CCLE_pseudo_bulk_file)
    #print(data.X)
    #print(data.X.shape)
    data_df = data.to_df()
    temp_df = data_df.T.merge(gene_length_df, left_index = True, right_index = True)    
    nm = norm()
    nm.tpm(df = temp_df, gl = 'gene_length')
    norm_temp_df = nm.tpm_norm
    ## TPM values (samples x genes)
    data_df_norm = norm_temp_df.T
    ## log2(TPM+1)
    data_df_log2 = np.log2(data_df_norm + 1)
    ## scaling
    data_mat = sparse.csr_matrix(data_df_log2.values)
    data.X = sample_scaling(data_mat)
    #print(data.X)
    data.write_h5ad(in_path + out_h5ad)
