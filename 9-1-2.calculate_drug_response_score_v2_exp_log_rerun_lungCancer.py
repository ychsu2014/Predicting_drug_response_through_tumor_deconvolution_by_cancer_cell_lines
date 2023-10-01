in_TCGA_pro_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/lung_cancer_all_oncoKB_19/"
in_TCGA_idx_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/idx_files/"
in_TCGA_drug_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/lung_cancer_all_oncoKB_19/"
in_PRISM = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/primary-screen-replicate-collapsed-logfold-change_filtered.txt"
in_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"

import os
import math

all_files = os.listdir(in_TCGA_pro_path)
pred_files = []
for afile in all_files:
    if "prediction.txt" in afile:
        pred_files.append(afile)

# for skin cancer rerun
#pred_files = ["SkinCancer_RNA_geneSubset_prediction.txt"]

# cell_line to depmap_ID
f_ID = open(in_ID)
lines = f_ID.readlines()
lines = lines[1:]
cell_ID_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    cell = cols[3]
    ID = cols[4]
    cell_ID_dict[cell] = ID
f_ID.close()

# key: ID, value: dict(key: drug_name, value: drug_response_value)
ID_drug_response_dict = {}
f = open(in_PRISM)
lines = f.readlines()
header = lines[0]
drug_list = header.strip("\n").split(",")
drug_list = drug_list[1:]
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split(",")
    ID = cols[0]
    drug_response = cols[1:]
    temp_dict = {}
    for i in range(len(drug_list)):
        temp_dict[drug_list[i]] = drug_response[i]
    ID_drug_response_dict[ID] = temp_dict
f.close()

for pred_file in pred_files:
    print(pred_file)
    f = open(in_TCGA_pro_path + pred_file)
    f_idx = open(in_TCGA_idx_path + pred_file.replace("_RNA_overlapped_genes_geneSubset_prediction.txt", "_RNA_idx_sID.txt"))
    #fout = open(in_TCGA_drug_path + pred_file.replace("_RNA_geneSubset_prediction.txt", "_drug_score.txt"), "w")
    fout2 = open(in_TCGA_drug_path + pred_file.replace("_RNA_overlapped_genes_geneSubset_prediction.txt", "_drug_score_None.txt"), "w")
    #fout.write("TCGA_sID\t" + "\t".join(drug_list) + "\n")
    fout2.write("TCGA_sID\t" + "\t".join(drug_list) + "\n")
    lines = f.readlines()
    header = lines[0]
    cell_line_list = header.strip("\n").split("\t")[1:]
    lines = lines[1:]
    # key: idx, value: TCGA sID
    idx_lines = f_idx.readlines()
    idx_sID_dict = {}
    for idx_line in idx_lines:
        cols = idx_line.strip("\n").split("\t")
        idx = cols[0]
        sID = cols[1]
        idx_sID_dict[idx] = sID
    # calculate drug response score
    for line in lines:
        cols = line.strip("\n").split("\t")
        idx = cols[0]
        TCGA_sID = idx_sID_dict[idx]
        #fout.write(TCGA_sID)
        fout2.write(TCGA_sID)
        pros = cols[1:]
        # Key: drug name, value: calculated score
        drug_score_dict = {}
        for adrug in drug_list:
            drug_score = 0
            NA_idx = False
            for i in range(len(pros)):
                cell = cell_line_list[i]
                proprotion = pros[i]
                CCLE_ID = cell_ID_dict[cell]
                drug_response_temp_dict = ID_drug_response_dict[CCLE_ID]
                drug_response = drug_response_temp_dict[adrug]
                if drug_response == "NA":
                    NA_idx = True
                else:
                    drug_score += float(proprotion) * (2 ** float(drug_response))
            #fout.write("\t" + str(drug_score))
            if NA_idx == True:
                fout2.write("\t" + str(None))
            else:
                fout2.write("\t" + str(math.log(drug_score, 2)))
        #fout.write("\n")
        fout2.write("\n")
    f.close()
    f_idx.close()
    #fout.close()
    fout2.close()
