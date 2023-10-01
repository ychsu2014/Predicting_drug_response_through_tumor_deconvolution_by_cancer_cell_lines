in_data = "/home/CBBI/hsuy1/projects/drugResponse/data/depmap_2022Q2_CCLE/CCLE_mutations.csv"
in_sID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/selected_cell_lines_for_lung_cancer_all_oncoKB_19.txt"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/lung_cancer_all_oncoKB_19/"

#in_sID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/selected_cell_lines_for_lung_cancer_oncoKB_with_drugs_11.txt"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/rerun_lung_cancer/lung_cancer_oncoKB_with_drugs_11/"

#in_sID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/selected_cell_lines_for_lung_cancer_all_oncoKB_mut_and_CGC_genes_30.txt"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/rerun_lung_cancer/lung_cancer_all_oncoKB_mut_and_CGC_genes_30/"

out_suffix = "_filtered_194.txt"
w_file_prefix = "LungCancer"

# sID_list
f = open(in_sID)
lines = f.readlines()
sID_list = []
for line in lines:
    sID = line.strip("\n")
    sID_list.append(sID)
f.close()

f = open(in_data)
lines = f.readlines()
header = lines[0]
lines = lines[1:]
fout = open(out_path + w_file_prefix + out_suffix, "w")
fout.write(header)
for line in lines:
    cols = line.strip("\n").split(",")
    sID = cols[15]
    if sID in sID_list:
        fout.write(line)
fout.close()
f.close()
