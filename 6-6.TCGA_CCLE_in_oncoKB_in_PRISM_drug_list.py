# make all drug name to lower case
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/"
in_PRISM = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_primary_treatment_info.txt"
in_TCGA = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7781.txt"
in_drug_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/"
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/in_PRISM_drug/TCGA_CCLE_oncoKB_PRISM_drug_oneDrugOneLine.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/in_PRISM_drug/TCGA_CCLE_oncoKB_PRISM_drug_count.txt"
out_file3 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/in_PRISM_drug/TCGA_CCLE_oncoKB_PRISM_drug_count_by_scRNA.txt"
out_file4 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/in_PRISM_drug/TCGA_CCLE_oncoKB_PRISM_drug_list.txt"

import os

# dict of PRISM drug information
f = open(in_PRISM)
lines = f.readlines()
lines = lines[1:]
PRISM_drug_info_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    drug_name = cols[2].lower()
    moa = cols[5]
    target = cols[6]
    disease_area = cols[7]
    indication = cols[8]
    PRISM_drug_info_dict[drug_name] = (moa, target, disease_area, indication)
f.close()

PRISM_drug_list = list(set(PRISM_drug_info_dict.keys()))
print(len(PRISM_drug_list))

# dict of cancer types
f = open(in_TCGA)
lines = f.readlines()
lines = lines[1:]
TCGA_cancer_type_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    sID = cols[0]
    TCGA_cancer_type = cols[1]
    scRNA_cancer_type = cols[2]
    TCGA_cancer_type_dict[sID] = (TCGA_cancer_type, scRNA_cancer_type)
f.close()

# get list of files containing TCGA and CCLE overlapped mutations that are in oncoKB
all_files = os.listdir(in_path)
oncoKB_files = []
for afile in all_files:
    if "_in_oncoKB_final.txt" in afile:
        oncoKB_files.append(afile)

# parse drug scores
all_files = os.listdir(in_drug_path)
drug_score_files = []
for afile in all_files:
    if "_addCancerType_addDrugName_None.txt" in afile:
        drug_score_files.append(afile)

TCGA_sID_drug_score = {}
for drug_score_file in drug_score_files:
    f = open(in_drug_path + drug_score_file)
    lines = f.readlines()
    drug_list = lines[1].strip("\n").split("\t")[3:]
    lines = lines[2:]
    for line in lines:
        cols = line.strip("\n").split("\t")
        sID = cols[0]
        drug_scores = cols[3:]
        for i in range(len(drug_list)):
            adrug = drug_list[i].lower()
            adrug_score = drug_scores[i]
            TCGA_sID_drug_score[(sID, adrug)] = adrug_score
    f.close()

# drug is in the last column
# one drug one line
# calculate number of drugs
fout1 = open(out_file1, "w")
fout2 = open(out_file2, "w")
fout3 = open(out_file3, "w")
fout4 = open(out_file4, "w")
# TCGA + oncoKB + drug + PRISM + CCLE cell lines
header_cols1 = ["TCGA_cancer_types", "scRNA_cancer_types", "sample", "chr", "start", "end", "reference", "alt", "gene", "effect", "Amino_Acid_Change", "DNA_VAF", "SIFT", "PolyPhen", "AA_change", "GRCh38_chr", "GRCh38_start", "GRCh38_end", "Level", "Gene", "Alterations", "Cancer Types", "oncoKB_drugs", "Drug", "PRISM_drug_moa", "PRISM_drug_target", "PRISM_drug_disease_area", "PRISM_drug_indication", "drug_score", "CCLE_cell_lines"]
header_cols2 = ["TCGA_cancer_types", "scRNA_cancer_types", "TCGA_sID", "drug_num", "drugs"]
fout1.write("\t".join(header_cols1) + "\n")
fout2.write("\t".join(header_cols2) + "\n")
fout4.write("TCGA_sID\tdrugs\tdrug_score\n")
TCGA_sID_drugs_dict = {}
TCGA_sID_PRISM_drugs_dict = {}
cancerType_drugs_dict = {}
cancerType_PRISM_drugs_dict = {}
for oncoKB_file in oncoKB_files:
    f = open(in_path + oncoKB_file)
    lines = f.readlines()
    for line in lines:
        cols = line.strip("\n").split("\t")
        sID = cols[0]
        drug_list = cols[len(cols)-1].split(", ")
        TCGA_line = cols[:16]
        cancer_types = TCGA_cancer_type_dict[sID]
        scRNA_cancer_type = cancer_types[1]
        oncoKB_line = cols[-5:]
        CCLE_cell_lines = list(set(cols) - set(TCGA_line) - set(oncoKB_line))
        for drug in drug_list:
            if drug != "":
                if drug.lower() in PRISM_drug_info_dict.keys():
                    PRISM_line = PRISM_drug_info_dict[drug.lower()]
                    wdrug_score = TCGA_sID_drug_score[(sID, drug.lower())]
                    fout1.write("\t".join(list(cancer_types)) + "\t" + "\t".join(TCGA_line) + "\t" + "\t".join(oncoKB_line) + "\t" + drug + "\t" + "\t".join(list(PRISM_line)) + "\t" + wdrug_score + "\t" + "\t".join(CCLE_cell_lines) + "\n")
                    if sID in TCGA_sID_PRISM_drugs_dict.keys():
                        old_list = TCGA_sID_PRISM_drugs_dict[sID]
                        old_list.append(drug)
                        old_list = list(set(old_list))
                        TCGA_sID_PRISM_drugs_dict[sID] = old_list
                    else:
                        TCGA_sID_PRISM_drugs_dict[sID] = [drug]
                    if scRNA_cancer_type in cancerType_PRISM_drugs_dict.keys():
                        old_list = cancerType_PRISM_drugs_dict[scRNA_cancer_type]
                        old_list.append(drug)
                        old_list = list(set(old_list))
                        cancerType_PRISM_drugs_dict[scRNA_cancer_type] = old_list
                    else:
                        cancerType_PRISM_drugs_dict[scRNA_cancer_type] = [drug]
                if sID in TCGA_sID_drugs_dict.keys():
                    old_list = TCGA_sID_drugs_dict[sID]
                    old_list.append(drug)
                    old_list = list(set(old_list))
                    TCGA_sID_drugs_dict[sID] = old_list
                else:
                    TCGA_sID_drugs_dict[sID] = [drug]
                if scRNA_cancer_type in cancerType_drugs_dict.keys():
                    old_list = cancerType_drugs_dict[scRNA_cancer_type]
                    old_list.append(drug)
                    old_list = list(set(old_list))
                    cancerType_drugs_dict[scRNA_cancer_type] = old_list
                else:
                    cancerType_drugs_dict[scRNA_cancer_type] = [drug]
    f.close()

for sID in TCGA_sID_drugs_dict.keys():
    drug_list = TCGA_sID_drugs_dict[sID]
    cancer_types = TCGA_cancer_type_dict[sID]
    if sID in TCGA_sID_PRISM_drugs_dict.keys():
        drug_list2 = TCGA_sID_PRISM_drugs_dict[sID]
    else:
        drug_list2 = []
    fout2.write("\t".join(cancer_types) + "\t" + sID + "\t" + str(len(drug_list)) + "\t" + "\t".join(drug_list) + "\t" + str(len(drug_list2)) + "\t" + "\t".join(drug_list2) + "\n")
    fout4.write(sID)
    for adrug in drug_list2:
        wdrug_score = TCGA_sID_drug_score[(sID, adrug.lower())]
        fout4.write("\t" + adrug + "\t" + wdrug_score)
    fout4.write("\n")

for scRNA_cancer_type in cancerType_drugs_dict.keys():
    drug_list_byCancer = cancerType_drugs_dict[scRNA_cancer_type]
    if scRNA_cancer_type in cancerType_PRISM_drugs_dict.keys():
        drug_list2_byCancer = cancerType_PRISM_drugs_dict[scRNA_cancer_type]
    else:
        drug_list2_byCancer = []
    fout3.write(scRNA_cancer_type + "\t" + str(len(drug_list_byCancer)) + "\t" + "\t".join(drug_list_byCancer) + "\t" + str(len(drug_list2_byCancer)) + "\t" + "\t".join(drug_list2_byCancer) + "\n")
    

fout1.close()
fout2.close()
fout3.close()
fout4.close()
