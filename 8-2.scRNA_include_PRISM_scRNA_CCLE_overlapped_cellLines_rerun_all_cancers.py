in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/greedy_algorithm_filtering_conditions_and_selected_cell_lines_for_all_cancers.txt"
in_ID_cell = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/preprocessing_onlyNormalized/by_cancertypes/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/preprocessing_onlyNormalized/"
cell_type_suffix = "_celltypes.txt"
count_suffix = "_norm_counts_all.txt"
include_conditions = ["oncoKB + CGC", "oncoKB", "oncoKB actionable target"]

import os

# subfolder dict: key: filtering conditions; value: subfolder
subfolder_dict = {}
subfolder_dict["oncoKB + CGC"] = "oncoKB_CGC"
subfolder_dict["oncoKB"] = "oncoKB"
subfolder_dict["oncoKB actionable target"] = "oncoKB_actionable_target"

# cell type files
all_files = os.listdir(in_path)
cell_type_files = []
for afile in all_files:
    if cell_type_suffix in afile:
        cell_type_files.append(afile)
cell_type_files.remove("Neuroblastoma_celltypes.txt")
cell_type_files.remove("Bone_Cancer_celltypes.txt")

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

# cell lines included
f = open(in_file)
lines = f.readlines()
lines = lines[1:]
# key: cancer_type (ex.ColonColorectalCancer), value: list of cell lines
cancer_cells_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    cancer_type = cols[0]
    cell_lines = cols[3]
    inc_condition = cols[1]
    exec("cancer_cells_dict[('" + cancer_type + "', '" + inc_condition + "')] = " + cell_lines)
f.close()

# filter the cell lines include in the cell type and count files
for cell_type_file in cell_type_files:
    cancer_type = cell_type_file.replace(cell_type_suffix, "").replace("_", "")
    count_file = cell_type_file.replace(cell_type_suffix, count_suffix)
    for inc_condition in include_conditions:
        out_subfolder = subfolder_dict[inc_condition]
        inc_cell_lines = cancer_cells_dict[(cancer_type, inc_condition)]
        if len(inc_cell_lines) >= 2:
            inc_cell_names = []
            for acell in inc_cell_lines:
                cell_name = ID_cell_dict[acell]
                inc_cell_names.append(cell_name)
            in_cell_type_file = in_path + cell_type_file
            in_count_file = in_path + count_file
            out_cell_type_file = out_path + out_subfolder + "/" + cell_type_file
            out_count_file = out_path + out_subfolder + "/" + count_file
            filterCellLines(in_cell_type_file, in_count_file, out_cell_type_file, out_count_file, inc_cell_names)
