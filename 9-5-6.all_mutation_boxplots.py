# draw boxplot for cancer-gene-drug combinations
in_mut_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/"
in_drug_info = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_primary_treatment_info.txt"
in_idx_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/idx_files/"
in_drug_via = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/drug_score_lung19/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/boxplots/all_mutations/"


idx_file_suffix = "_RNA_idx_sID.txt"
via_file_suffix = "_addCancerType_addDrugName_None.txt"

import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import os

#filter_cancer_type = "BladderCancer"
#filter_gene = "RP1"
#filter_aa_change = "R677*"
#filter_drug = "LY2606368"
#filter_p_value = "5.45E-113"

#filter_cancer_type = "EndometrialUterineCancer"
#filter_gene = "ANKRD55"
#filter_aa_change = "R327Q"
#filter_drug = "A-674563"
#filter_p_value = "1.01E-145"

#filter_cancer_type = "ColonColorectalCancer"
#filter_gene = "ZNF331"
#filter_aa_change = "M1?"
#filter_drug = "KW-2478"
#filter_p_value = "4.02E-59"

#filter_cancer_type = "ColonColorectalCancer"
#filter_gene = "MYPOP"
#filter_aa_change = "P380Hfs*26"
#filter_drug = "KW-2478"
#filter_p_value = "1.65E-58"

#filter_cancer_type = "ColonColorectalCancer"
#filter_gene = "MEGF8"
#filter_aa_change = "T1603Qfs*60"
#filter_drug = "KW-2478"
#filter_p_value = "1.84E-72"

#filter_cancer_type = "ColonColorectalCancer"
#filter_gene = "IGFN1"
#filter_aa_change = "R3593Gfs*8"
#filter_drug = "KW-2478"
#filter_p_value = "1.84E-72"

#filter_cancer_type = "HeadandNeckCancer"
#filter_gene = "FAT1"
#filter_aa_change = "S3373*"
#filter_drug = "evodiamine"
#filter_p_value = "3.90E-187"

#filter_cancer_type = "HeadandNeckCancer"
#filter_gene = "FAT1"
#filter_aa_change = "S3373*"
#filter_drug = "talazoparib"
#filter_p_value = "5.23E-157"

#filter_drug = "lestaurtinib"
#filter_p_value = "1.68E-179"

#filter_cancer_type = "BreastCancer"
#filter_gene = "PI16"
#filter_aa_change = "V284I"
#filter_drug = "chloropyramine"
#filter_p_value = "9.68E-190"

#filter_drug = "ONX-0914"
#filter_p_value = "6.74E-176"

#filter_drug = "OTS167"
#filter_p_value = "6.58E-300"

#filter_drug = "puromycin"
#filter_p_value = "2.91E-181"

filter_cancer_type = "BreastCancer"
filter_gene = "PIK3CA"
filter_aa_change = "H1047R"
filter_drug = "alpelisib"
filter_p_value = "0.0114"

filter_aa_change = "E545K"
filter_p_value = "0.00208"

def getFiles(inPath, inStr, inCancer):
    allFiles = os.listdir(inPath)
    outFiles = []
    for afile in allFiles:
        if inStr in afile and inCancer in afile:
            outFiles.append(afile)
    return(outFiles)

def drawBoxplot(inData1, inData2, inLabel1, inLabel2, inTitle, inPvalue, outFile):
    ##### modify the position if needed
    text_x = 0.75
    text_y = 0.4
    #####
    all_fontsize = 14
    data = [inData1, inData2]
    plt.figure(figsize =(5, 5))
    # Creating axes instance
    #ax = fig.add_axes([0.2, 0.2, 0.8, 0.8])
    # Creating plot
    plt.boxplot(data)
    #plt.scatter([1] * len(inData1), inData1, facecolors='none', edgecolors='black')
    data1Num = str(len(inData1))
    data2Num = str(len(inData2))
    #plt.set_xticklabels(["CR (n = " + CRnum + ")", "SD (n = " + SDnum + ")"], fontsize = 16)
    #ax.yaxis.set_tick_params(labelsize = 16)
    plt.xticks([1,2], [inLabel1 + "(n = " + data1Num + ")", inLabel2 + "(n = " + data2Num + ")"], fontsize = all_fontsize, rotation = 45)
    plt.yticks(fontsize = all_fontsize)
    plt.title(inTitle, fontsize = all_fontsize)
    plt.ylabel("Predicted viability", fontsize = all_fontsize)
    plt.text(text_x, text_y, "adjusted p-value: " + str(inPvalue))
    plt.tight_layout()
    plt.savefig(outFile, format = "pdf")
    #plt.show()
    plt.close()

# key: cancer type, value: list of sIDs
all_ID_files = getFiles(in_idx_path, idx_file_suffix, filter_cancer_type)
ID_files = []
for ID_file in all_ID_files:
    if "_not_in_all_RNA_idx_sID.txt" not in ID_file and "_RNA_idx_sID.txt" in ID_file:
        ID_files.append(ID_file)

#cancer_sID_dict = {}
cancer_all_sID_list = []
for ID_file in ID_files:
    f = open(in_idx_path + ID_file)
    lines = f.readlines()
    for line in lines:
        cols = line.strip("\n").split("\t")
        TCGA_sID = cols[1]
        cancer_all_sID_list.append(TCGA_sID)
    f.close()
print("Number of samples: " + str(len(cancer_all_sID_list)))
print("ID completed.")

# predicted cell viability
# key: (sID, drug), value: predicted_via (float, excluded None)
via_files = getFiles(in_drug_via, via_file_suffix, filter_cancer_type)
sID_drug_via_dict = {}
all_drug_list = []
for via_file in via_files:
    print(via_file)
    f = open(in_drug_via + via_file)
    lines = f.readlines()
    header = lines[1]
    drug_list = header.strip("\n").split("\t")[3:]
    lines = lines[2:]
    for line in lines:
        cols = line.strip("\n").split("\t")
        sID = cols[0]
        predicted_vias = cols[3:]
        for i in range(len(drug_list)):
            predicted_via = predicted_vias[i]
            drug = drug_list[i]
            # BRD-K28042756-001-01-9::5::HTS
            if drug == "":
                continue
            if predicted_via != "None":
                all_drug_list.append(drug)
                sID_drug_via_dict[(sID, drug)] = float(predicted_via)
    f.close()

all_drug_list = list(set(all_drug_list))
print("check if blank in all_drug_list:")
print("" in all_drug_list)
print("Number of drugs:" + str(len(all_drug_list)))
print("Viability completed.")

# TCGA mutation sID dict
# key: mutation, value: list of TCGA sID
# mutation: oncoKB: (gene, proteinChange); CGC: gene
mut_file = filter_cancer_type + "_GRCh38_7781filtered.txt"
mut_df = pd.read_csv(in_mut_path + mut_file, sep = "\t")
gene_idx = set(np.where(mut_df["gene"] == filter_gene)[0])
aa_change_idx = set(np.where(mut_df["AA_change"] == filter_aa_change)[0])
select_idx = list(gene_idx.intersection(aa_change_idx))
mut_sID_list = list(mut_df.iloc[select_idx]["sample"])
cancer_mut_sID_dict = {}
cancer_mut_sID_dict[filter_gene + "_" + filter_aa_change] = mut_sID_list

print("Mutation completed.")

cancer_type_dict = {}
cancer_type_dict["BileDuctCancer"] = "Bile Duct Cancer"
cancer_type_dict["BladderCancer"] = "Bladder Cancer"
cancer_type_dict["BrainCancer"] = "Brain Cancer"
cancer_type_dict["BreastCancer"] = "Breast Cancer"
cancer_type_dict["ColonColorectalCancer"] = "Colon/Colorectal Cancer"
cancer_type_dict["EndometrialUterineCancer"] = "Endometrial/Uterine Cancer"
cancer_type_dict["EsophagealCancer"] = "Esophageal Cancer"
cancer_type_dict["GastricCancer"] = "Gastric Cancer"
cancer_type_dict["HeadandNeckCancer"] = "Head and Neck Cancer"
cancer_type_dict["KidneyCancer"] = "Kidney Cancer"
cancer_type_dict["LiverCancer"] = "Liver Cancer"
cancer_type_dict["LungCancer"] = "Lung Cancer"
cancer_type_dict["OvarianCancer"] = "Ovarian Cancer"
cancer_type_dict["PancreaticCancer"] = "Pancreatic Cancer"
cancer_type_dict["ProstateCancer"] = "Prostate Cancer"
cancer_type_dict["Sarcoma"] = "Sarcoma"
cancer_type_dict["SkinCancer"] = "Skin Cancer"
cancer_type_dict["ThyroidCancer"] = "Thyroid Cancer"

print(len(cancer_mut_sID_dict.keys()))
count = 0
for agene in cancer_mut_sID_dict.keys():
    count += 1
    if count % 100 == 0:
        print(count)
    mut_sID_list = cancer_mut_sID_dict[agene]
    control_sID_list = list(set(cancer_all_sID_list).difference(set(mut_sID_list)))
    if len(mut_sID_list) > 1 and len(control_sID_list) > 1:
        ###### modified for drawing boxplot
        mut_value_list = []
        con_value_list = []
        for mut_sID in mut_sID_list:
            if (mut_sID, filter_drug) in sID_drug_via_dict.keys():
                mut_via = sID_drug_via_dict[(mut_sID, filter_drug)]
                mut_value_list.append(mut_via)
        for con_sID in control_sID_list:
            if (con_sID, filter_drug) in sID_drug_via_dict.keys():
                con_via = sID_drug_via_dict[(con_sID, filter_drug)]
                con_value_list.append(con_via)
        fig_cancer_type = cancer_type_dict[filter_cancer_type]
        out_file = out_path + filter_cancer_type + "_" + filter_gene + "_" + filter_aa_change + "_" + filter_drug + ".pdf"
        drawBoxplot(mut_value_list, con_value_list, agene + "-mut", "wildtype", filter_drug + " in " + fig_cancer_type, filter_p_value, out_file)
