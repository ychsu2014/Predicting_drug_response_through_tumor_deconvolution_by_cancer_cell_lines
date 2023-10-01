in_celltype = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/pan_cancer_celltypes.txt"
in_exp = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/pan_cancer_norm_counts_all.txt"
in_cancertype = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/ID_cellline_cancertype_v2_manuallyChecked_noHeader.txt"
in_IDnum = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/ID_num.txt"
out_path = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/by_cancertypes_all_cell_lines"

import pandas as pd

# ID_num, Cell_line
df_celltype = pd.read_csv(in_celltype, sep = "\t")
map_dict = {}
map_dict[df_celltype.columns[0]] = "ID_num"
map_dict[df_celltype.columns[1]] = "Cell_line"
df_celltype.rename(columns = map_dict, inplace = True)

# ID, Cell_line, Cancer_type
df_cacnertype = pd.read_csv(in_cancertype, sep = "\t", header = None, names = ["ID", "Cell_line", "Cancer_type"])

# ID_num, ID
df_ID = pd.read_csv(in_IDnum, sep = "\t", names = ["ID_num", "ID"])

# merge
df_merge = df_celltype.merge(df_ID, on = "ID_num", how = "left")
df_merge = df_merge.merge(df_cacnertype, on = "ID", how = "left")

# gene expression data
df_exp = pd.read_csv(in_exp, sep = "\t")
exp_cols = ["ID_num"] + list(df_exp.columns)[1:]
map_dict2 = {}
for i in range(len(df_exp.columns)):
    map_dict2[df_exp.columns[i]] = exp_cols[i]
df_exp.rename(columns = map_dict2, inplace = True)

df_merge = df_merge.merge(df_exp, on = "ID_num", how = "left")

# Index(['ID_num', 'Cell_line_x', 'ID', 'Cell_line_y', 'Cancer_type', 'A1BG', 'A2M', 'AARD', 'ABCA5', ...'ECSCR', 'PIK3R3', 'RGS5'], dtype='object', length=2224)

# split by cancer type
cancer_type_list = list(set(df_merge["Cancer_type"]))
for cancer_type in cancer_type_list:
    print(cancer_type)
    out_prefix = cancer_type.replace(" ", "_").replace("/", "_")
    # cell types (cell lines)
    cell_type_data = df_merge[df_merge["Cancer_type"] == cancer_type]
    cell_type_data = cell_type_data[["ID_num", "Cell_line_x"]]
    cell_type_data = cell_type_data.rename(columns = {"ID_num": "", "Cell_line_x": "Celltype"})
    cell_type_data.to_csv(out_path + "/" + out_prefix + "_celltypes.txt", sep = "\t", index = None)
    df_merge[df_merge["Cancer_type"] == cancer_type].iloc[:,5:].to_csv(out_path + "/" + out_prefix + "_norm_counts_all.txt", sep = "\t", index = True, header = True)
