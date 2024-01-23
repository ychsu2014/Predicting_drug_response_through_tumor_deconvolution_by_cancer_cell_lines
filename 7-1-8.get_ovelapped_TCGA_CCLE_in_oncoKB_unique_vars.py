in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/"
#in_ID_cell = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
out_file = in_path + "TCGA_CCLE_in_oncoKB_uni_var.txt"
out_file2 = in_path + "TCGA_CCLE_in_oncoKB_uni_var_stat.txt"
out_file3 = in_path + "TCGA_CCLE_in_oncoKB_uni_entry_with_duplicated_vars.txt"
out_file4 = in_path + "TCGA_CCLE_in_oncoKB_uni_entry_with_duplicated_vars_v2.txt"

import os

# depmap_ID to cell line name
#f = open(in_ID_cell)
#lines = f.readlines()
#lines = lines[1:]
#ID_cell_dict = {}
#for line in lines:
#    cols = line.strip("\n").split("\t")
#    ID = cols[4]
#    cell_line = cols[3]
#    ID_cell_dict[ID] = cell_line
#f.close()

# files for TCGA/CCLE overlapped mutations that are also in oncoKb
all_files = os.listdir(in_path)
in_files = []
for afile in all_files:
    if "_in_oncoKB_final.txt" in afile:
        in_files.append(afile)

fout = open(out_file, "w")
fout2 = open(out_file2, "w")
fout3 = open(out_file3, "w")
fout4 = open(out_file4, "w")
for in_file in in_files:
    f = open(in_path + in_file)
    lines = f.readlines()
    uni_var_list = []
    # list of unique entries (may contain duplicated vars)
    uni_entry_list = []
    uni_entry_list2 = []
    for line in lines:
        cols = line.strip("\n").split("\t")
        var = "\t".join(cols[1:len(cols)-5]) + "\n"
        uni_var_list.append(var)
        #entry = "\t".join(cols[1:]) + "\n"
        entry2 = "\t".join(cols[:16]) + "\t" + "\t".join(cols[len(cols)-5:]) + "\t" + "\t".join(cols[16:len(cols)-5]) + "\n"
        uni_entry_list.append(line)
        uni_entry_list2.append(entry2)
    uni_var_list = list(set(uni_var_list))
    uni_entry_list = list(set(uni_entry_list))
    uni_entry_list2 = list(set(uni_entry_list2))
    fout2.write(in_file + "\t" + str(len(uni_var_list)) + "\n")
    for uni_var in uni_var_list:
        fout.write(in_file + "\t" + uni_var)
    for uni_entry in uni_entry_list:
        fout3.write(in_file + "\t" + uni_entry)
    for uni_entry in uni_entry_list2:
        fout4.write(in_file + "\t" + uni_entry)
fout.close()
fout2.close()
fout3.close()
fout4.close()
