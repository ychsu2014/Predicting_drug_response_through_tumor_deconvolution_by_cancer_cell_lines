in_file = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/all_RNA.h5ad"
in_ID = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/all_RNA_idx_sID.txt"
in_phe = "/home/yuching/projects/drugResponse/data/TCGA/TCGA_scRNA_cancer_type_mapping_updated20230228.txt"
in_phe_ID = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/TCGA_pheno_RNA_DNA_overlapped_sID_cancerType_list.txt"
out_file = "/home/yuching/projects/drugResponse/data/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7781.txt"
out_path1 = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/by_cancertypes/"
out_path2 = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/idx_files/"
out1_suffix = "_RNA.h5ad"
out2_suffix = "_RNA_idx_sID.txt"
out3_suffix = "_sID_not_in_all_RNA_idx_sID.txt"


import anndata as ad

def parseToDict(inFile, keyColNum, valueColNum):
    f = open(inFile)
    lines = f.readlines()
    lines = lines[1:]
    outDict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        keyCol = cols[keyColNum]
        valueCol = cols[valueColNum]
        if keyCol != "":
            if keyCol in outDict.keys():
                old_list = outDict[keyCol]
                old_list.append(valueCol)
                old_list = list(set(old_list))
                outDict[keyCol] = old_list
            else:
                outDict[keyCol] = [valueCol]
    f.close()
    return(outDict)

scRNACanerType_TCGACancerType_dict = parseToDict(in_phe, 2, 0)
#update_only_skin_cancer = {}
#update_only_skin_cancer["Skin Cancer"] = scRNACanerType_TCGACancerType_dict["Skin Cancer"]
TCGACancerType_sID_dict = parseToDict(in_phe_ID, 3, 0)

f_ID = open(in_ID)
lines = f_ID.readlines()
sID_idx_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    idx = cols[0]
    sID = cols[1]
    sID_idx_dict[sID] = idx
f_ID.close()

fout = open(out_file, "w")
fout.write("TCGA_ID\tTCGA_cancer_type\tscRNA_cancer_type\n")
scCancerType_idx_dict = {}
for sc_cancer_type in scRNACanerType_TCGACancerType_dict.keys(): #update_only_skin_cancer.keys():
    print(sc_cancer_type)
    TCGA_cancer_types = scRNACanerType_TCGACancerType_dict[sc_cancer_type] #update_only_skin_cancer[sc_cancer_type]
    keep_idx_list = []
    write_scRNA_cancer_type = sc_cancer_type
    if " " in write_scRNA_cancer_type:
        write_scRNA_cancer_type = write_scRNA_cancer_type.replace(" ", "")
    if "/" in write_scRNA_cancer_type:
        write_scRNA_cancer_type = write_scRNA_cancer_type.replace("/", "")
    fout2 = open(out_path2 + write_scRNA_cancer_type + out2_suffix, "w")
    fout3 = open(out_path2 + write_scRNA_cancer_type + out3_suffix, "w")
    for TCGA_cancer_type in TCGA_cancer_types:
        print(TCGA_cancer_type)
        sID_list = TCGACancerType_sID_dict[TCGA_cancer_type]
        for sID in sID_list:
            fout.write(sID + "\t" + TCGA_cancer_type + "\t" + sc_cancer_type + "\n")
            if sID in sID_idx_dict.keys():
                idx = sID_idx_dict[sID]
                keep_idx_list.append(idx)
                keep_idx_list = list(keep_idx_list)
                fout2.write(idx + "\t" + sID + "\n")
            else:
                fout3.write(sID + "\n")
    scCancerType_idx_dict[sc_cancer_type] = keep_idx_list
fout.close()
fout2.close()
fout3.close()


all_RNA = ad.read_h5ad(in_file)
for cancer_type in scCancerType_idx_dict.keys():
    print(cancer_type)
    idx_list = scCancerType_idx_dict[cancer_type]
    w_scRNA_cancer_type = cancer_type
    if " " in w_scRNA_cancer_type:
        w_scRNA_cancer_type = w_scRNA_cancer_type.replace(" " , "")
    if "/" in w_scRNA_cancer_type:
        w_scRNA_cancer_type = w_scRNA_cancer_type.replace("/", "")
    RNA_temp = all_RNA[idx_list]
    RNA_temp.write_h5ad(out_path1 + w_scRNA_cancer_type + out1_suffix)
