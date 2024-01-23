in_clinical_path = "/home/yuching/projects/drugResponse/data/TCGA/clinical_data/"
in_predict_path = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/drug_score/exponential_log/drug_score_lung19_collapsed/"
in_mapping = "/home/yuching/projects/drugResponse/data/TCGA/TCGA_scRNA_cancer_type_mapping_updated20230228_addedTCGAabbreviation.txt"
in_drug_mapping = "/home/yuching/projects/drugResponse/data/TCGA/TCGA_clinical_drug_mapped_to_PRISM_drugname.txt"
out_path = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/clinical/res_nonRes_boxplots_v3/"
out_p = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/clinical/res_nonRes_boxplots_v3/p_values.txt"

drug_response = ["Complete Response", "Partial Response"]
drug_non_response = ["Stable Disease", "Clinical Progressive Disease"]

import os
import matplotlib.pyplot as plt
import rpy2
import rpy2.robjects as robjects
import numpy as np
import math
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# dict for TCGA cancer type and abbreviation mapping
# key: TCGA cancer abbreviation, value: TCGA cancer type
# ex. CHOL, cholangiocarcinoma
f = open(in_mapping)
lines = f.readlines()
lines = lines[1:]
TCGA_abbrev_cancerType_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    cancer_type = cols[0]
    cancer_abbrev = cols[1]
    TCGA_abbrev_cancerType_dict[cancer_abbrev] = cancer_type
f.close()

# dict for predicted viability
# key: (TCGA_cancer_type, drug), value: a dict (key: TCGA_sID, value: predicted_via)
all_files = os.listdir(in_predict_path)
predict_files = []
for afile in all_files:
    if "_addCancerType_addDrugName_None_collapsed.txt" in afile:
        predict_files.append(afile)

##### for lung cancer boxplot
#predict_files = ["LungCancer_addCancerType_addDrugName_None_collapsed.txt"]   

cancer_drug_sID_score_dict = {}
for predict_file in predict_files:
    f = open(in_predict_path + predict_file)
    lines = f.readlines()
    header = lines[1]
    drug_list = header.strip("\n").split("\t")[3:]
    lines = lines[2:]
    for line in lines:
        cols = line.strip("\n").split("\t")
        TCGA_sID = cols[0].replace("-01", "")
        TCGA_cancer_type = cols[1]
        drug_scores = cols[3:]
        for i in range(len(drug_scores)):
            drug_name = drug_list[i]#.lower()
            drug_score = drug_scores[i]
            if (TCGA_cancer_type, drug_name) in cancer_drug_sID_score_dict.keys():
                old_dict = cancer_drug_sID_score_dict[(TCGA_cancer_type, drug_name)]
                old_dict[TCGA_sID] = drug_score
                cancer_drug_sID_score_dict[(TCGA_cancer_type, drug_name)] = old_dict
            else:
                temp_dict = {}
                temp_dict[TCGA_sID] = drug_score
                cancer_drug_sID_score_dict[(TCGA_cancer_type, drug_name)] = temp_dict
    f.close()

# TCGA drug name to PRISM drug name
# TCGA_PRISM_drug_dict: key: TCGA drug name, value: list of PRISM drug names
# 109 unique PRISM drugs
f = open(in_drug_mapping)
TCGA_PRISM_drug_dict = {}
lines = f.readlines()
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    TCGA_drug_name_temp = cols[0]
    PRISM_drug_name_temp = cols[1]
    if TCGA_drug_name_temp in TCGA_PRISM_drug_dict.keys():
        old_list = TCGA_PRISM_drug_dict[TCGA_drug_name_temp]
        old_list.append(PRISM_drug_name_temp)
        old_list = list(set(old_list))
        TCGA_PRISM_drug_dict[TCGA_drug_name_temp] = old_list
    else:
        TCGA_PRISM_drug_dict[TCGA_drug_name_temp] = [PRISM_drug_name_temp]
f.close()

# key: (TCGA_cancer_abbrev, drug), value: a dict (key: "Complete Response" or "Stable Disease", value: [(TCGA_sID, predicted_viability)]
#clinical_files = os.listdir(in_clinical_path)
clinical_files = ["LUSC_clinical_drug.txt", "ESCA_clinical_drug.txt"]
cancer_drug_response_sID_score_dict = {}
for clinical_file in clinical_files:
    f = open(in_clinical_path + clinical_file)
    TCGA_cancer_abbrev = clinical_file.strip("_clinical_drug.txt")
    TCGA_cancer_type = TCGA_abbrev_cancerType_dict[TCGA_cancer_abbrev]
    lines = f.readlines()
    lines = lines[1:]
    for line in lines:
        cols = line.replace('"', '').strip("\n").split("\t")
        response = cols[20]
        if response in drug_response or response in drug_non_response:
            drug_name = cols[14]
            if drug_name in TCGA_PRISM_drug_dict.keys():
                PRISM_drug_name_list = TCGA_PRISM_drug_dict[drug_name]
                TCGA_sID = cols[1]
                for PRISM_drug_name in PRISM_drug_name_list:
                    if (TCGA_cancer_type, PRISM_drug_name) in cancer_drug_sID_score_dict.keys():
                        prediction_dict = cancer_drug_sID_score_dict[(TCGA_cancer_type, PRISM_drug_name)]
                        if TCGA_sID in prediction_dict.keys():
                            predict_score = prediction_dict[TCGA_sID]
                            if predict_score != "None":
                                if (TCGA_cancer_abbrev, PRISM_drug_name) in cancer_drug_response_sID_score_dict.keys():
                                    old_dict = cancer_drug_response_sID_score_dict[(TCGA_cancer_abbrev, PRISM_drug_name)]
                                    if response in old_dict.keys():
                                        old_list = old_dict[response]
                                        old_list.append((TCGA_sID, predict_score))
                                        old_list = list(set(old_list))
                                        old_dict[response] = old_list
                                    else:
                                        old_dict[response] = [(TCGA_sID, predict_score)]
                                    cancer_drug_response_sID_score_dict[(TCGA_cancer_abbrev, PRISM_drug_name)] = old_dict
                                else:
                                    temp_dict = {}
                                    temp_dict[response] = [(TCGA_sID, predict_score)]
                                    cancer_drug_response_sID_score_dict[(TCGA_cancer_abbrev, PRISM_drug_name)] = temp_dict
    f.close()

print("complete.")
# boxplot + p-value (using R package)
# inTitle: "tamoxifen in BRCA (p-value = 0.01)"
#def drawBoxplot(inCRdata, inSDdata, inTitle, outFile):
#    data = [inCRdata, inSDdata]
#    fig = plt.figure(figsize =(5, 3))
#    # Creating axes instance
#    ax = fig.add_axes([0, 0, 1, 1])
#    # Creating plot
#    bp = ax.boxplot(data)
#    CRnum = str(len(inCRdata))
#    SDnum = str(len(inSDdata))
#    ax.set_xticklabels(["CR (n = " + CRnum + ")", "SD (n = " + SDnum + ")"], fontsize = 16)
#    ax.yaxis.set_tick_params(labelsize = 16)
#    plt.title(inTitle)
#    plt.ylabel("predicted values", fontsize = 16)
#    plt.savefig(outFile)

def drawBoxplot(inCRdata, inSDdata, inTitle, outFile, markDots = False):
    data = [inCRdata, inSDdata]
    fig = plt.figure(figsize =(10, 7))
    # Creating axes instance
    #ax = fig.add_axes([0.2, 0.2, 0.8, 0.8])
    # Creating plot
    plt.boxplot(data)
    CRnum = str(len(inCRdata))
    SDnum = str(len(inSDdata))
    ### mark all dots when the number of samples in each group is small (n < 10)
    if markDots == True:
        x_list = [1] * int(CRnum) + [2] * int(SDnum)
        y_list = inCRdata + inSDdata
        plt.scatter(x_list, y_list)
        outFile = outFile.replace(".pdf", "_v2.pdf")
    #plt.set_xticklabels(["CR (n = " + CRnum + ")", "SD (n = " + SDnum + ")"], fontsize = 16)
    #ax.yaxis.set_tick_params(labelsize = 16)
    plt.xticks([1,2], ["Responder (n = " + CRnum + ")", "Non-responder (n = " + SDnum + ")"], fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.title(inTitle, fontsize = 20)
    plt.ylabel("Predicted values", fontsize = 20)
    plt.tight_layout()
    #plt.show()
    plt.savefig(outFile, format = "pdf")
    plt.close()

def ttestByR(inCRdata, inSDdata):
    CRvector = robjects.FloatVector(inCRdata)
    SDvector = robjects.FloatVector(inSDdata)
    t_test = robjects.r["t.test"]
    result = t_test(CRvector, SDvector)
    p_value = float(np.array(result[2])[0])
    return(p_value)

# key: (TCGA_cancer_abbrev, drug), value: a dict (key: "Complete Response" or "Stable Disease", value: [(TCGA_sID, predicted_viability)]
fout = open(out_p, "w")
for (TCGA_cancer_abbrev, drug) in cancer_drug_response_sID_score_dict.keys():
    response_dict = cancer_drug_response_sID_score_dict[(TCGA_cancer_abbrev, drug)]
    if ("Complete Response" in response_dict.keys() or "Partial Response" in response_dict.keys()) and ("Stable Disease" in response_dict.keys() or "Clinical Progressive Disease" in response_dict.keys()):
        CR_list = []
        SD_list = []
        if "Complete Response" in response_dict.keys():
            CR_list += response_dict["Complete Response"]
        if "Partial Response" in response_dict.keys():
            CR_list += response_dict["Partial Response"]
        if "Stable Disease" in response_dict.keys():
            SD_list += response_dict["Stable Disease"]
        if "Clinical Progressive Disease" in response_dict.keys():
            SD_list += response_dict["Clinical Progressive Disease"]
        if len(CR_list) > 1 and len(SD_list) > 1:
            print(TCGA_cancer_abbrev)
            print(drug)
            CR_value_list = []
            SD_value_list = []
            for aCR in CR_list:
                predict_via = aCR[1]
                CR_value_list.append(float(predict_via))
            for aSD in SD_list:
                predict_via = aSD[1]
                SD_value_list.append(float(predict_via))
            print(CR_value_list)
            print(SD_value_list)
            p_value = round(ttestByR(CR_value_list, SD_value_list), 4)
            fout.write(drug + "\t" + TCGA_cancer_abbrev + "\t" + str(p_value) + "\n")
            in_title = drug + " in " + TCGA_cancer_abbrev + " (p-value = " + str(p_value) + ")" 
            out_file = out_path + drug + "_" + TCGA_cancer_abbrev + "_" + str(p_value) + ".pdf"
            drawBoxplot(CR_value_list, SD_value_list, in_title, out_file, True)
fout.close()
