# include only cell lines included in scRNA, CCLE, PRISM primary dataset
#in_BRCA_cellLines = ["HCC1428", "HDQP1", "EFM192A", "HCC38", "MDAMB436", "ZR751", "BT474", "HMC18", "HCC1419", "BT549", "MCF7", "T47D", "CAMA1"]
#in_cell_lines = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/scRNA_cell_lines.txt"
cell_line_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/"
in_lung_cell_line_files = ["selected_cell_lines_for_lung_cancer_all_oncoKB_19.txt", "selected_cell_lines_for_lung_cancer_oncoKB_with_drugs_11.txt", "selected_cell_lines_for_lung_cancer_all_oncoKB_mut_and_CGC_genes_30.txt"]
in_ID_cell = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/by_cancertypes_selected_cell_lines/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/"
cell_type_suffix = "_celltypes.txt"
count_suffix = "_norm_counts_all.txt"
cancer_type_prefix = "Lung_Cancer"

import os

# create folders for the three conditions
new_folder_list = []
for cell_file in in_lung_cell_line_files:
    new_folder = cell_file.replace("selected_cell_lines_for_", "").replace(".txt", "")
    new_folder_list.append(new_folder)
    os.mkdir(out_path + new_folder)

# depmap ID to cell line names
f = open(in_ID_cell)
lines = f.readlines()
lines = lines[1:]
ID_cell_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    cell_line = cols[3]
    depmap_ID = cols[4]
    ID_cell_dict[depmap_ID] = cell_line
f.close()

def filterCellLines(inCelltypes, inCounts, outCelltypes, outCounts, includeCells):
    # cell lines data by cancer type
    fout_cell_types = open(outCelltypes, "w")
    f_cell_types = open(inCelltypes)
    lines = f_cell_types.readlines()
    header = lines[0]
    fout_cell_types.write(header)
    lines = lines[1:]
    included_idx_list = []
    for line in lines:
        cols = line.strip("\n").split("\t")
        idx = cols[0]
        cell_line = cols[1].split("_")[0]
        if cell_line in includeCells:
            fout_cell_types.write(idx + "\t" + cell_line + "\n")
            included_idx_list.append(idx)
    fout_cell_types.close()
    f_cell_types.close()
    # count data by cancer type
    fout_counts = open(outCounts, "w")
    f_counts = open(inCounts)
    lines = f_counts.readlines()
    header = lines[0]
    fout_counts.write(header)
    lines = lines[1:]
    for line in lines:
        cols = line.strip("\n").split("\t")
        idx = cols[0]
        if idx in included_idx_list:
            fout_counts.write(line)
    fout_counts.close()
    f_counts.close()

for cell_file in in_lung_cell_line_files:
    sub_folder = cell_file.replace("selected_cell_lines_for_", "").replace(".txt", "")
    in_cell_type_file = in_path + cancer_type_prefix + cell_type_suffix
    in_count_file = in_path + cancer_type_prefix + count_suffix
    out_cell_type_file = out_path + sub_folder + "/" + cancer_type_prefix + cell_type_suffix
    out_count_file = out_path + sub_folder + "/" + cancer_type_prefix + count_suffix
    f = open(cell_line_path + cell_file)
    lines = f.readlines()
    include_cell_line_list = []
    for line in lines:
        cell_ID = line.strip("\n")
        cell_name = ID_cell_dict[cell_ID]
        include_cell_line_list.append(cell_name)
    print(out_cell_type_file)
    print(include_cell_line_list)
    filterCellLines(in_cell_type_file, in_count_file, out_cell_type_file,out_count_file,include_cell_line_list)
    f.close()
