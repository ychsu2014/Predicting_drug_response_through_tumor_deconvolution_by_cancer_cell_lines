# input the file path to simulation data & the filename of the simulation data

import scanpy as sc
import pandas as pd
import numpy as np
import sys

data_path = sys.argv[1]
simu_name = sys.argv[2]
out_file = sys.argv[3]

in_true_file = data_path + simu_name[:-5] + "_processed.h5ad"
in_pred_file = data_path + simu_name[:-5] + "_prediction.txt"

fout = open(out_file, "a")

true_h5ad = sc.read_h5ad(in_true_file)
t_cell_types = true_h5ad.uns["cell_types"]

pred_df = pd.read_csv(in_pred_file, sep = "\t", index_col = 0).T
p_cell_types = pred_df.index

if set(t_cell_types) != set(p_cell_types):
    print("true cell types:")
    print(set(t_cell_types))
    print("prediction cell types:")
    print(set(p_cell_types))
    sys.exit("Cell types inconsistent between true proportions and predicted proportions. Please check.")

all_cell_types = list(set(t_cell_types))
sIDs = pred_df.columns

true_pred_data_dict = {}

def appendToDict(inValue, inType, inDict = true_pred_data_dict):
    if inType in inDict.keys():
        oldList = inDict[inType]
        oldList.append(inValue)
        inDict[inType] = oldList
    else:
        inDict[inType] = [inValue]
    return(inDict)

print(len(sIDs))
count = 0
for sID in sIDs:
    count += 1
    if count % 1000 == 0:
        print(count)
    for i in range(len(all_cell_types)):
        cell_type = all_cell_types[i]
        # for author's simulated data
        #temp_true = true_h5ad.obs[cell_type].loc[sID]
        # for my own preprocessed data
        #temp_true = true_h5ad.obs[cell_type].iloc[sID]
        # for pan cancer data
        temp_true = true_h5ad.obs[cell_type].loc[str(sID)]
        temp_pred = pred_df.loc[cell_type, sID]
        true_pred_data_dict = appendToDict(temp_true, "true_y")
        true_pred_data_dict = appendToDict(temp_pred, "pred_y")

#all_keys = list(true_pred_data_dict.keys())
#print(true_pred_data_dict[all_keys[0]])


true_pred_df = pd.DataFrame.from_dict(true_pred_data_dict)
#true_pred_df.to_csv(out_file, header=True, index=None, sep='\t')

y_true = true_pred_df["true_y"]
y_pred = true_pred_df["pred_y"]
cor = np.corrcoef(y_true, y_pred)[0][1]
mean_true = np.mean(y_true)
mean_pred = np.mean(y_pred)
var_true = np.var(y_true)
var_pred = np.var(y_pred)
sd_true = np.std(y_true)
sd_pred = np.std(y_pred)
numerator = 2 * cor * sd_true * sd_pred
denominator = var_true + var_pred + (mean_true - mean_pred) ** 2
ccc = numerator / denominator
fout.write(simu_name + "\t" + str(ccc) + "\n")
print(simu_name)
print(ccc)