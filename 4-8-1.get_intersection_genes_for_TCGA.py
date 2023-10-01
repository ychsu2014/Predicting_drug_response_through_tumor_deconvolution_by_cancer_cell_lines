# get intersection genes and scaling to [0,1]
# Since the TCGA RNA data is already log2 transformed, we can not use sceden process which includes log2 transformation.

# for TCGA data
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"
### for testing the correctness of this script
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/test_script/"
# for CCLE data
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"
in_train_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"

import os
import anndata as ad
from scipy import sparse
from sklearn import preprocessing as pp

def sample_scaling(x):
    mms = pp.MinMaxScaler(feature_range=(0, 1), copy=True)
    # it scales features so transpose is needed
    x = mms.fit_transform(x.T.toarray()).T
    x = sparse.csr_matrix(x)
    return(x)

train_files = os.listdir(in_train_path)
train_processed_h5ad_files = []
for train_file in train_files:
    if "training" in train_file and "processed.h5ad" in train_file:
        train_processed_h5ad_files.append(train_file)

### for testing the correctness of this script 
#train_processed_h5ad_files = ["training_Bile_Duct_Cancer_processed.h5ad"]

# for rerun skin cancer
#rerun_files = ["training_Skin_Cancer_processed.h5ad"]

for train_file in train_processed_h5ad_files: #rerun_files:
    in_file = train_file.replace("training", "").replace("processed.h5ad", "").replace("_", "") + "_RNA_overlapped_genes.h5ad"
    in_file = in_path + in_file
    in_train = in_train_path + train_file
    out_h5ad = in_file.replace("_RNA_overlapped_genes.h5ad", "_RNA_overlapped_genes_geneSubset.h5ad")
    adata = ad.read_h5ad(in_file)
    train = ad.read_h5ad(in_train)
    adata_gene_list = list(adata.var_names)
    train_gene_list = list(train.var_names)
    overlapped_gene_list = list(set(adata_gene_list).intersection(set(train_gene_list)))
    # in the Scaden preprocess command, they do the gene subsetting first, then the scaling
    # include overlapped genes between training dataset and TCGA data
    adata_sub = adata[:,overlapped_gene_list]
    # scaling
    adata_sub.X = sample_scaling(adata_sub.X)
    adata_sub.write_h5ad(out_h5ad)
