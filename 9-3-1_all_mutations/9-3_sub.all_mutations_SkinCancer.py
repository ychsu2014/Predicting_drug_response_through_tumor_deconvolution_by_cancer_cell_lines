###### modified for parallele processing
out_cancer_type = "SkinCancer"

in_mut = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/" + out_cancer_type + "_GRCh38_7781filtered.txt"
in_drug_info = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_primary_treatment_info.txt"

in_idx_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/idx_files/"
in_drug_via = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/"
###### modified for parallele processing
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/mut_grouped_drug_prediction_t_test_p_values_all_mutations_" + out_cancer_type + "_v4.txt"

idx_file_suffix = "_RNA_idx_sID.txt"
via_file_suffix = "_addCancerType_addDrugName_None.txt"

import os
import rpy2.robjects as robjects
import numpy as np
import pandas as pd

###### modified for parallele processing
def getFiles(inPath, inStr, inCancer):
    allFiles = os.listdir(inPath)
    outFiles = []
    for afile in allFiles:
        if inStr in afile and inCancer in afile:
            outFiles.append(afile)
    return(outFiles)

# key: cancer type, value: list of sIDs
all_ID_files = getFiles(in_idx_path, idx_file_suffix, out_cancer_type)
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
via_files = getFiles(in_drug_via, via_file_suffix, out_cancer_type)
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
f = open(in_mut)
lines = f.readlines()
lines = lines[1:]
cancer_mut_sID_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    AA_change = cols[12]
    gene = cols[6]
    sID = cols[0]
    mut_comb = (gene, AA_change)
    if mut_comb in cancer_mut_sID_dict.keys():
        old_list = cancer_mut_sID_dict[mut_comb]
        old_list.append(sID)
        old_list = list(set(old_list))
        cancer_mut_sID_dict[mut_comb] = old_list
    else:
        cancer_mut_sID_dict[mut_comb] = [sID]
f.close()

print("Number of mutation: " + str(len(cancer_mut_sID_dict.keys())))
print("Mutation completed.")

# drug target dict (for on target prediction)
# key: drug, value: list of genes
#f = open(in_drug_info)
#lines = f.readlines()
#drug_target_dict = {}
#for line in lines:
#    cols = line.strip("\n").split("\t")
#    drug_name = cols[2]
#    targets = cols[6].split(", ")
#    drug_target_dict[drug_name] = targets
#f.close()

def addToDict(inDict, inKey, in2ndKey, inValue):
    if inKey in inDict.keys():
        tempDict = inDict[inKey]
        if in2ndKey in tempDict.keys():
            dataList = tempDict[in2ndKey]
            dataList.append(inValue)
            dataList = list(set(dataList))
            tempDict[in2ndKey] = dataList
        else:
            tempDict[in2ndKey] = [inValue]
        inDict[inKey] = tempDict
    else:
        tempDict = {}
        tempDict[in2ndKey] = [inValue]
        inDict[inKey] = tempDict
    return(inDict)

# all data dict
# key: (cancer_type, gene, drug), value: dict(key: 0 (controls), or 1 (cases); value: [(predicted_via)])
#drug_via_file = in_drug_via + via_files[0]
#print(drug_via_file)
#drug_df = pd.read_csv(drug_via_file, sep = "\t", low_memory = False)
#old_header = list(drug_df.columns)
#new_header = list(drug_df.iloc[0])
#rename_dict = {}
#for i in range(len(new_header)):
#    rename_dict[old_header[i]] = new_header[i]
#drug_df.rename(columns = rename_dict, inplace = True)
#drug_df.drop(0, inplace = True)
#drug_df.set_index("TCGA_sID", inplace = True)

print(len(cancer_mut_sID_dict.keys()))
all_data_dict = {}
count = 0
for mut_comb in cancer_mut_sID_dict.keys():
    count += 1
    if count % 100 == 0:
        print(count)
    mut_sID_list = cancer_mut_sID_dict[mut_comb]
    control_sID_list = list(set(cancer_all_sID_list).difference(set(mut_sID_list)))
    if len(mut_sID_list) > 1 and len(control_sID_list) > 1:
        for drug in all_drug_list:
            mut_value_list = []
            con_value_list = []
            for mut_sID in mut_sID_list:
                if (mut_sID, drug) in sID_drug_via_dict.keys():
                    mut_via = sID_drug_via_dict[(mut_sID, drug)]
                    mut_value_list.append(mut_via)
            for con_sID in control_sID_list:
                if (con_sID, drug) in sID_drug_via_dict.keys():
                    con_via = sID_drug_via_dict[(con_sID, drug)]
                    con_value_list.append(con_via)
            #mut_value_list = drug_df.loc[mut_sID_list][drug]
            #con_value_list = drug_df.loc[control_sID_list][drug]
            temp_dict = {}
            temp_dict[1] = mut_value_list
            temp_dict[0] = con_value_list
            all_data_dict[(mut_comb, drug)] = temp_dict

print("All data completed.")

def ttestByR(inData1, inData2):
    vector1 = robjects.FloatVector(inData1)
    vector2 = robjects.FloatVector(inData2)
    t_test = robjects.r["t.test"]
    result = t_test(vector1, vector2)
    p_value = float(np.array(result[2])[0])
    return(p_value)

# key: (cancer_type, (gene, p_hgvs), drug), value: dict(key: 0 (controls), or 1 (cases); value: [(predicted_via)])
# dataframe with:
# (1) differences of predicted log2 fold change between mutant and wildtype
# (2) p-values
# (3) adjusted p-values 
# (4) cancer types
# (5) mutation
# (6) drug

fout = open(out_file, "w")
fout.write("Cancer_type\tgene\tprotein_hgvs\tdrug\tmut_num\twt_num\tmut_mean\twt_mean\tmut_wt_diff\tp_values\n")
print(len(all_data_dict.keys()))
count = 0
for (mut_comb, drug) in all_data_dict.keys():
    (gene, aa_change) = mut_comb
    count += 1
    if count % 10000 == 0:
        print(count)
    grouped_via_dict = all_data_dict[(mut_comb, drug)]
    if 1 in grouped_via_dict.keys() and 0 in grouped_via_dict.keys():
        mut_via_list = grouped_via_dict[1]
        wt_via_list = grouped_via_dict[0]
        if len(mut_via_list) > 1 and len(wt_via_list) > 1:
            mut_num = len(mut_via_list)
            wt_num = len(wt_via_list)
            mut_mean = np.mean(mut_via_list)
            wt_mean = np.mean(wt_via_list)
            mut_wt_diff = mut_mean - wt_mean
            p_value = ttestByR(mut_via_list, wt_via_list)
            fout.write("\t".join([out_cancer_type, gene, aa_change, drug, str(mut_num), str(wt_num), str(mut_mean), str(wt_mean), str(mut_wt_diff), str(p_value)]) + "\n")
fout.close()
