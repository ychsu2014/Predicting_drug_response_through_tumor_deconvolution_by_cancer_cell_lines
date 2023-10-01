# get the cancer type from metadata
in_UMI = "/home/yuching/projects/drugResponse/data/scRNA/raw_data_SCP542/UMIcount_data.txt"
in_meta = "/home/yuching/projects/drugResponse/data/scRNA/raw_data_SCP542/Metadata.txt"
out_file = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/ID_cellline_cancertype.txt"

# UMI
f = open(in_UMI)
lines = f.readlines()
ID_list = lines[0].strip("\n").split("\t")[1:]
cell_line_list = lines[1].strip("\n").split("\t")[1:]
print(len(ID_list))
print(len(cell_line_list))
ID_cellline_dict_UMI = {}
for i in range(len(ID_list)):
    ID = ID_list[i]
    cell_line = cell_line_list[i]
    ID_cellline_dict_UMI[ID] = cell_line
print(i)
f.close()

# metadata
f = open(in_meta)
lines = f.readlines()
lines = lines[2:]
ID_cellline_dict_meta = {}
ID_cancertype_dict_meta = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    ID = cols[0]
    cell_line = cols[1]
    cancer_type = cols[3]
    ID_cellline_dict_meta[ID] = cell_line
    ID_cancertype_dict_meta[ID] = cancer_type
    if cell_line == "HS729_SOFT_TISSUE":
        print(cancer_type)
f.close()

# check all the ID in metadata is also in the UMI header
not_in_num = 0
for ID_meta in ID_cellline_dict_meta.keys():
    if ID_meta not in ID_cellline_dict_UMI.keys():
        print(ID_meta)
        not_in_num += 1
print("Num of IDs not in UMI header: " + str(not_in_num))

fout = open(out_file, "w")
fout.write("ID\tcell_line_UMI\tcell_line_meta\tcancer_type_meta\n")
for ID in ID_cellline_dict_UMI.keys():
    cell_line_UMI = ID_cellline_dict_UMI[ID]
    if ID in ID_cellline_dict_meta.keys():
        cell_line_meta = ID_cellline_dict_meta[ID]
        cancer_type_meta = ID_cancertype_dict_meta[ID]
    else:
        cell_line_meta = ""
        cancer_type_meta = ""
    fout.write("\t".join([ID, cell_line_UMI, cell_line_meta, cancer_type_meta]) + "\n")
fout.close()
