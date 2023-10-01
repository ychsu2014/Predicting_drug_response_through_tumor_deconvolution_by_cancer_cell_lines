# TCGA_ID TCGA_cancer_type        scRNA_cancer_type
#TCGA-W5-AA2Z-01 cholangiocarcinoma      Bile Duct Cancer
in_cancer = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7781.txt"
in_drug_name = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_primary_treatment_info.txt"
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/"
#in_suffix = "_drug_score.txt"
# drug matrix added cancer types and drug names
#out_file1_suffix = "_addCancerType_addDrugName.txt"
# TCGA_sID, TCGA_cancer_type, scRNA_cancer_type, drugs(passed threshold)
#out_file2_suffix = "_drugsPassedThreshold.txt"
in_suffix = "_drug_score_None.txt"
out_file1_suffix = "_addCancerType_addDrugName_None.txt"
out_file2_suffix = "_drugsPassedThreshold_None.txt"
out_file = in_path + "TCGA_sID_drugCount_passed_threshold_0.3.txt"

in_threshold = 0.3

import os
import numpy as np
import pandas as pd
import ast

# key: TCGA sID, value: (TCGA_cancer_type, scRNA_cancer_type)
f = open(in_cancer)
lines = f.readlines()
TCGA_sID_cancerType_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    sID = cols[0]
    TCGA_cancer = cols[1]
    scRNA_cancer = cols[2]
    TCGA_sID_cancerType_dict[sID] = (TCGA_cancer, scRNA_cancer)
f.close()

# key: drug ID, value: drug name
f = open(in_drug_name)
lines = f.readlines()
PRISM_drug_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    drug_ID = cols[0]
    drug_name = cols[2]
    PRISM_drug_dict[drug_ID] = drug_name
f.close()

# drug score files
all_files = os.listdir(in_path)
drug_score_files = []
for afile in all_files:
    if in_suffix in afile:
        drug_score_files.append(afile)

# write to file
fout = open(out_file, "w")
fout.write("TCGA_sID\tTCGA_cancer_type\tscRNA_cancer_type\tdrug_counts\n")
for drug_file in drug_score_files:
    f = open(in_path + drug_file)
    print(drug_file)
    fout1 = open(in_path + drug_file.replace(in_suffix, out_file1_suffix), "w")
    fout2 = open(in_path + drug_file.replace(in_suffix, out_file2_suffix), "w")
    lines = f.readlines()
    # write the header lines for output file 1
    header = lines[0]
    header_cols = header.strip("\n").split("\t")
    new_header_cols = [header_cols[0]] + ["TCGA_cancer_type", "scRNA_cancer_type"] + header_cols[1:]
    fout1.write("\t".join(new_header_cols) + "\n")
    new_header_cols2 = [header_cols[0]] + ["TCGA_cancer_type", "scRNA_cancer_type"]
    for drug_ID in header_cols[1:]:
        drug_name = PRISM_drug_dict[drug_ID]
        new_header_cols2.append(drug_name)
    fout1.write("\t".join(new_header_cols2) + "\n")
    drug_list = new_header_cols2[3:]
    drug_series = pd.Series(drug_list)
    # write the header for output file 2
    fout2.write("\t".join(["TCGA_sID", "TCGA_cancer_type", "scRNA_cancer_type", "number_of_drugs", "drug_names"]) + "\n")
    lines = lines[1:]
    for line in lines:
        cols = line.strip("\n").split("\t")
        sID = cols[0]
        drug_scores = cols[1:]
        # write drug scores
        (TCGA_cancer, scRNA_cancer) = TCGA_sID_cancerType_dict[sID]
        fout1.write(sID + "\t" + TCGA_cancer + "\t" + scRNA_cancer + "\t" + "\t".join(drug_scores) + "\n")
        # write drugs passed threshold
        drug_score_float = list(map(ast.literal_eval, drug_scores))
        drug_score_series = pd.Series(drug_score_float)
        passed_drug_list = list(drug_series.loc[np.where(drug_score_series < in_threshold)])
        fout2.write(sID + "\t" + TCGA_cancer + "\t" + scRNA_cancer + "\t" + str(len(passed_drug_list)) + "\t" + "\t".join(passed_drug_list) + "\n")
        fout.write(sID + "\t" + TCGA_cancer + "\t" + scRNA_cancer + "\t" + str(len(passed_drug_list)) + "\n")
    f.close()
    fout1.close()
    fout2.close()
fout.close()
