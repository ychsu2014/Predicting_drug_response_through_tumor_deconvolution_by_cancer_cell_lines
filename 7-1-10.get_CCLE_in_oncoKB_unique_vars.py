in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/"
in_file = in_path + "CCLE_all_mutations_filtered194_in_oncoKB_final.txt"
out_file = in_path + "CCLE_all_mutations_filtered194_in_oncoKB_uni_var.txt"
out_file2 = in_path + "CCLE_all_mutations_filtered194_in_oncoKB_uni_var_stat.txt"
out_file3 = in_path + "CCLE_all_mutations_filtered194_in_oncoKB_uni_entry_with_duplicated_vars.txt"

import os

fout = open(out_file, "w")
fout2 = open(out_file2, "w")
fout3 = open(out_file3, "w")

def parseToDict(inDict, inCancerType, inData):
    if inCancerType in inDict.keys():
        dataList = inDict[inCancerType]
        dataList.append(inData)
        dataList = list(set(dataList))
        inDict[inCancerType] = dataList
    else:
        inDict[inCancerType] = [inData]
    return(inDict)

f = open(in_file)
lines = f.readlines()
uni_var_dict = {}
# list of unique entries (may contain duplicated vars)
uni_entry_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    scRNA_cancer_type = cols[32]
    var_cols = cols[:15] + cols[16:len(cols)-5]
    #var = "\t".join(cols[1:len(cols)-5]) + "\n"
    var = "\t".join(var_cols) + "\n"
    uni_var_dict = parseToDict(uni_var_dict, scRNA_cancer_type, var)
    uni_entry_dict = parseToDict(uni_entry_dict, scRNA_cancer_type, line)

for cancer_type in uni_var_dict.keys():
    var_list = uni_var_dict[cancer_type]
    fout2.write(cancer_type + "\t" + str(len(var_list)) + "\n")
    for uni_var in var_list:
        fout.write(cancer_type + "\t" + uni_var)

for cancer_type in uni_entry_dict.keys():
    entry_list = uni_entry_dict[cancer_type]
    for uni_entry in entry_list:
        fout3.write(cancer_type + "\t" + uni_entry)

fout.close()
fout2.close()
fout3.close()
