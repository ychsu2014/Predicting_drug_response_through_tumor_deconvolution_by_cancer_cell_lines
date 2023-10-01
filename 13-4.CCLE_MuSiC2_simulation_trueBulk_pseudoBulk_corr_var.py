true_folder = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"
# Scaden simulation
#pred_folder = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/"
# MuSiC2 simulation
pred_folder = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/h5ad/"
in_ID_cell = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/CCLE_all_RNA_idx_depmapID.txt"
out_file = pred_folder + "CCLE_cell_line_RNA_trueBulk_pseudoBulk_corrCoef_unexpVar_MuSiC2.txt"

num_of_genes_included = 10000

import os
import anndata as ad
from scipy import stats
import pandas as pd
#import numpy as np

all_files = os.listdir(pred_folder)
pred_files = []
for afile in all_files:
    if "_pseudo_bulk_norm_scaled.h5ad" in afile:
        pred_files.append(afile)

#pred_files = ["Thyroid_Cancer_pseudo_bulk_norm_scaled.h5ad"]

# key: idx, value: cell line
f = open(in_ID_cell)
lines = f.readlines()
lines = lines[1:]
idx_cell_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    idx = cols[0]
    cell_name = cols[1]
    idx_cell_dict[idx] = cell_name
f.close()

def square(x):
    return(x ** 2)

fout = open(out_file, "w")
fout.write("\t".join(["cancer_type", "sID", "cell_line", "correlation_coefficient", "p_value", "unexplained_variance"]) + "\n")
for pred_file in pred_files:
    cancer_type = pred_file.replace("_pseudo_bulk_norm_scaled.h5ad", "")
    true_file = cancer_type.replace("_","") + "_RNA_ovelapped_genes_geneSubset.h5ad"
    pred_data = ad.read_h5ad(pred_folder + pred_file)
    true_data = ad.read_h5ad(true_folder + true_file)
    # reset the index of pred_data (ex. pred_data: 0,1,2,3; true_data: "777", "857", "1117", "315")
    #pred_data.obs.set_index("sample_ID", inplace = True)
    all_sIDs = list(true_data.obs.index)
    true_df = true_data.to_df()
    pred_df = pred_data.to_df()
    for sID in all_sIDs:
        cell_name = idx_cell_dict[sID]
        ### using all genes
        #true_temp = list(true_df.loc[sID,])
        #pred_temp = list(pred_df.loc[sID, true_df.columns])
        ### using highly expressed genes
        true_temp_df = pd.DataFrame(true_df.loc[sID,])
        pred_temp_df = pred_df.loc[cell_name, true_df.columns]
        merge_df = true_temp_df.merge(pred_temp_df, left_index = True, right_index = True)
        merge_df.rename(columns = {merge_df.columns[0]:"true", merge_df.columns[1]: "pred"}, inplace = True)
        merge_df = merge_df[:num_of_genes_included]
        true_temp = list(merge_df["true"])
        pred_temp = list(merge_df["pred"])
        res = stats.pearsonr(true_temp, pred_temp)
        corr_coef = res[0]
        p_value = res[1]
        unexp_var = 1 - (res[0] ** 2)
        fout.write("\t".join([cancer_type, sID, cell_name, str(corr_coef), str(p_value), str(unexp_var)]) + "\n")
        #true_pred_diff = np.array(true_temp) - np.array(pred_temp)
        #diff_square = list(map(square, true_pred_diff))
        #unexp_var = sum(diff_square)
fout.close()
