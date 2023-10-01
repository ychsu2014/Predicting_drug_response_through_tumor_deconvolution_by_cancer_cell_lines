in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/pseudo_bulk_RNA/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/pseudo_bulk_RNA/quantile_normalization/"

import os
import anndata as ad
from sklearn import preprocessing as pp
from scipy import sparse
import qnorm
import numpy as np

all_files = os.listdir(in_path)
pseudo_bulk_files = []
for afile in all_files:
    if "_pseudo_bulk.h5ad" in afile:
        pseudo_bulk_files.append(afile)

pseudo_bulk_files = ["Bile_Duct_Cancer_pseudo_bulk.h5ad"]

def sample_scaling(x):
    mms = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
    ## it scales features so transpose is needed
    x = mms.fit_transform(x.T.toarray()).T
    x = sparse.csr_matrix(x)
    return(x)

for pseudo_bulk_file in pseudo_bulk_files:
    print(pseudo_bulk_file)
    out_h5ad = pseudo_bulk_file.replace(".h5ad", "_norm_scaled.h5ad")
    data = ad.read_h5ad(in_path + pseudo_bulk_file)
    data_df = data.to_df().T
    # quantile normalization (input: genes x samples)
    norm_data_df = qnorm.quantile_normalize(data_df, axis = 1)
    # log2(norm_value + 1)
    log2_norm_df = np.log2(norm_data_df + 1)
    # scaling (input: samples x genes)
    data_mat = sparse.csr_matrix(log2_norm_df.T.values)
    data.X = sample_scaling(data_mat)
    data.write_h5ad(out_path + out_h5ad)
