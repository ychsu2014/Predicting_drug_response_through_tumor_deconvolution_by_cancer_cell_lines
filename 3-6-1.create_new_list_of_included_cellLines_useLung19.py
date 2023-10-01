in_all_cell = "/home/yuching/projects/drugResponse/data/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
in_lung = "/home/yuching/projects/drugResponse/data/scRNA/selected_cell_lines_for_lung_cancer_all_oncoKB_19.txt"
out_file = "/home/yuching/projects/drugResponse/data/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list_v2_lung19.txt"

filter_cancer = "Lung Cancer"

# get the selected 19 lung cancer cell lines
f = open(in_lung)
lines = f.readlines()
selected_lung_cell_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    selected_lung_cell_list.append(cols[0])
f.close()

# filter out the unselected lung cancer cell lines, and write to new file
f = open(in_all_cell)
fout = open(out_file, "w")
lines = f.readlines()
header = lines[0]
fout.write(header)
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    ID = cols[4]
    cancer = cols[1]
    if cancer == filter_cancer:
        if ID in selected_lung_cell_list:
            fout.write(line)
    else:
        fout.write(line)
f.close()
fout.close()
