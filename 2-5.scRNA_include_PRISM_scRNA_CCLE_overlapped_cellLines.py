# keep only cell lines included in scRNA, CCLE, PRISM primary dataset
#in_BRCA_cellLines = ["HCC1428", "HDQP1", "EFM192A", "HCC38", "MDAMB436", "ZR751", "BT474", "HMC18", "HCC1419", "BT549", "MCF7", "T47D", "CAMA1"]
in_cell_lines = "/home/yuching/projects/drugResponse/data/scRNA/scRNA_cell_lines.txt"
in_path = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/by_cancertypes_all_cell_lines/"
out_path = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/by_cancertypes_selected_cell_lines/"
cell_type_suffix = "_celltypes.txt"
count_suffix = "_norm_counts_all.txt"

import os

all_files = os.listdir(in_path)

# cancertype_cell_dict: key: cancer types, value: list of cell lines
f_cell = open(in_cell_lines)
lines = f_cell.readlines()
header = lines[0]
lines = lines[1:]
cancertype_cell_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    cancer_type = cols[1]
    cell_line = cols[2]
    if cancer_type in cancertype_cell_dict.keys():
        old_list = cancertype_cell_dict[cancer_type]
        old_list.append(cell_line)
        old_list = list(set(old_list))
        cancertype_cell_dict[cancer_type] = old_list
    else:
        cancertype_cell_dict[cancer_type] = [cell_line]
f_cell.close()

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

for cancer_type in cancertype_cell_dict.keys():
    cancer_type_prefix = cancer_type
    if " " in cancer_type:
        cancer_type_prefix = cancer_type_prefix.replace(" ", "_")
    if "/" in cancer_type:
        cancer_type_prefix = cancer_type_prefix.replace("/", "_")
    temp = cancer_type_prefix + cell_type_suffix
    if temp in all_files:
        in_cell_type_file = in_path + cancer_type_prefix + cell_type_suffix
        in_count_file = in_path + cancer_type_prefix + count_suffix
        out_cell_type_file = out_path + cancer_type_prefix + cell_type_suffix
        out_count_file = out_path + cancer_type_prefix + count_suffix
        include_cell_line_list = cancertype_cell_dict[cancer_type]
        print(out_cell_type_file)
        print(include_cell_line_list)
        filterCellLines(in_cell_type_file, in_count_file, out_cell_type_file,out_count_file,include_cell_line_list)
