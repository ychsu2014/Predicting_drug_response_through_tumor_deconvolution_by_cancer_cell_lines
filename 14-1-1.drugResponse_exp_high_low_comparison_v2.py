in_CGC_file = "/home/UTHSCSA_drug/data/COSMIC/Census_allMon_Jan_30_21_08_13_2023.tsv"
in_RNA_path = "/home/UTHSCSA_drug/data/TCGA/RNA/by_cancertypes/"
in_RNA_lung_path = "/home/UTHSCSA_drug/data/TCGA/RNA/lung_cancer_all_oncoKB_19/"
in_drug_path = "/home/UTHSCSA_drug/data/TCGA/drug_score_lung19/"
in_idx_path = "/home/UTHSCSA_drug/data/TCGA/idx_files/"
out_file = "/home/UTHSCSA_drug/data/TCGA/exp_lung19/exp_high_low_drug_score_t_test_pvalue_" + in_cancer_type + "_part" + str(part_num) +".txt"
out_log = "/home/UTHSCSA_drug/data/TCGA/exp_lung19/exp_high_low_drug_score_t_test_pvalue_" + in_cancer_type + "_log_part" + str(part_num) +".txt"
out_error = "/home/UTHSCSA_drug/data/TCGA/exp_lung19/exp_high_low_drug_score_t_test_pvalue_" + in_cancer_type + "_error_part" + str(part_num) +".txt"
out_folder = "/data5/yuching"

drug_files = [in_cancer_type + "_addCancerType_addDrugName_None.txt"]

fout = open(out_file, "w")
fout_log = open(out_log, "w")
fout_error = open(out_error, "w")

# ThyroidCancer_RNA_overlapped_genes_geneSubset.h5ad
RNA_suffix = "_RNA_overlapped_genes_geneSubset.h5ad"
drug_suffix = "_addCancerType_addDrugName_None.txt"
idx_suffix = "_RNA_idx_sID.txt"

import os
import anndata as ad
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from pydash import at
from datetime import datetime
import math

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

def ttestByR(inData1, inData2):
    vector1 = robjects.FloatVector(inData1)
    vector2 = robjects.FloatVector(inData2)
    #t_test = robjects.r["t.test"]
    #result = t_test(vector1, vector2)
    stats = importr("stats")
    result = stats.t_test(vector1, vector2)
    p_value = float(np.array(result[2])[0])
    return(p_value)

# CGC list
f = open(in_CGC_file)
lines = f.readlines()
lines = lines[1:]
CGC_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    CGC_list.append(gene)
CGC_list = CGC_list[CGC_start_num:CGC_end_num]
fout_log.write("CGC parsed. " + str(datetime.now()) + "\n")

#drug_files = os.listdir(in_drug_path)
### for test
#drug_files = ["BileDuctCancer_addCancerType_addDrugName_None.txt"]
# header of the output file
fout.write("Cancer_type\tGene\tDrug_ID\tDrug_name\tHigh_exp_num\tLow_exp_num\tHigh_TCGA_ID\tLow_TCGA_ID\tHigh_group_DR\tLow_group_DR\tHigh_mean_DR\tLow_mean_DR\tHigh_low_DR_diff\tP_value\n")
for drug_file in drug_files:
    print(drug_file)
    cancer_type = drug_file.replace(drug_suffix, "")
    cancer_type = cancer_type_dict[cancer_type]
    # idx_ID_dict: key: idx, value: TCGA ID
    idx_file = drug_file.replace(drug_suffix, idx_suffix)
    f = open(in_idx_path + idx_file)
    lines = f.readlines()
    idx_ID_dict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        idx_ID_dict[cols[0]] = cols[1]
    print("idx file completed.")
    fout_log.write("idx file completed. " + str(datetime.now()) + "\n")
    # RNA_dict: key: idx, value: expression value
    RNA_file = drug_file.replace(drug_suffix, RNA_suffix)
    if "LungCancer" in RNA_file:
        RNA_data = ad.read_h5ad(in_RNA_lung_path + RNA_file)
    else:
        RNA_data = ad.read_h5ad(in_RNA_path + RNA_file)
    RNA_df = RNA_data.to_df()
    RNA_dict = RNA_df.to_dict()
    print("RNA file completed.")
    fout_log.write("RNA file completed. " + str(datetime.now()) + "\n")
    # drugID_drugName_dict: key: drug ID, value: drug name
    # drug_dict: key: drug, value:dicts: key: TCGA_ID, drug response
    drug_df = pd.read_csv(in_drug_path + drug_file, sep = "\t", low_memory = False)
    #drug_df.columns = drug_df.iloc[1]
    drug_df.set_index("TCGA_sID", inplace = True)
    drugID_drugName_df = drug_df[:1]
    drugID_drugName_df = drugID_drugName_df.iloc[:,2:]
    # BRD-M63173034-001-03-6::2.64076472::MTS004  clonixin-lysinate
    drugID_drugName_dict = drugID_drugName_df.to_dict("records")[0]
    drug_df.drop(["TCGA_sID"], inplace = True)
    drug_df.drop(["TCGA_cancer_type", "scRNA_cancer_type"], axis = 1, inplace = True)
    drug_dict = drug_df.to_dict()
    print("drug file completed.")
    fout_log.write("drug file completed. " + str(datetime.now()) + "\n")
    # for each gene, get the median and split the samples into high/low groups
    gene_list = list(RNA_dict.keys())
    ### for test
    #gene_list = ["A1BG"]
    count = 0
    for agene in CGC_list:
        if agene in gene_list:
            fout_log.write("start " + str(agene) + " " + str(datetime.now()) + "\n")
            count += 1
            if count % 10 == 0:
                fout_log.write(str(count) + "\t" + str(datetime.now()) + "\n")
            oneGene_RNA_dict = RNA_dict[agene]
            exp_values = list(oneGene_RNA_dict.values())
            exp_median = np.median(exp_values)
            high_group_ID_list = []
            low_group_ID_list = []
            fout_log.write("Start splitting TCGA IDs to high/low groups. " + str(datetime.now()) + "\n")
            # split to high/low groups
            for idx in oneGene_RNA_dict.keys():
                exp_value = oneGene_RNA_dict[idx]
                TCGA_ID = idx_ID_dict[idx]
                if exp_value <= exp_median:
                    low_group_ID_list.append(TCGA_ID)
                else:
                    high_group_ID_list.append(TCGA_ID)
            fout_log.write("Splitting ends. " + str(datetime.now()) + "\n")
            # get drug response from high/low groups
            high_group_str = "','".join(high_group_ID_list)
            low_group_str = "','".join(low_group_ID_list)
            drug_list = list(drug_dict.keys())
            ### for test
            fout_log.write("Start parsing drug response data. " + str(datetime.now()) + "\n")
            #drug_list = ["BRD-A00055058-001-01-0::2.325889319::MTS004"]
            #drug_list = ["BRD-A00077618-236-07-6::2.5::HTS"]
            drug_count = 0
            for adrug in drug_list:
                drug_count += 1
                if drug_count % 500 == 0:
                    print(str(count) + "-" + str(drug_count))
                drug_name = drugID_drugName_dict[adrug]
                fout_log.write("Start splitting drug response data. " + str(datetime.now()) + "\n")
                exec("high_group_drug = at(drug_dict['" + adrug + "'],'" + high_group_str + "')")
                exec("low_group_drug = at(drug_dict['" + adrug + "'],'" + low_group_str + "')")
                if "None" in high_group_drug or "None" in low_group_drug or np.nan in high_group_drug or np.nan in low_group_drug or "nan" in high_group_drug or "nan" in low_group_drug:
                    continue
                elif type(high_group_drug[0]) == float:
                    if np.isnan(high_group_drug[0]) or np.isnan(low_group_drug[0]) or math.isnan(high_group_drug[0]) or math.isnan(low_group_drug[0]):
                        continue
                    else:
                        continue
                        fout_error.write("\t".join([in_cancer_type, agene, str(high_group_drug), str(low_group_drug)]) + "\n")
                else:
                    high_group_drug = list(map(float, high_group_drug))
                    low_group_drug = list(map(float, low_group_drug))
                    high_drug_mean = np.mean(high_group_drug)
                    low_drug_mean = np.mean(low_group_drug)
                    high_low_drug_diff = high_drug_mean - low_drug_mean
                    fout_log.write("Start t-test. " + str(datetime.now()) + "\n")
                    p_value = ttestByR(high_group_drug, low_group_drug)
                    #fout.write("Cancer_type\tGene\tDrug_ID\tDrug_name\tHigh_exp_num\tLow_exp_num\tHigh_TCGA_ID\tLow_TCGA_ID\tHigh_mean_DR\tLow_mean_DR\tHigh_low_DR_diff\tP_value\n")
                    fout_log.write("Start writing data. " + str(datetime.now()) + "\n")
                    fout.write("\t".join([cancer_type, str(agene), str(adrug), str(drug_name), str(len(high_group_drug)), str(len(low_group_drug)), str(high_group_ID_list), str(low_group_ID_list), str(high_group_drug), str(low_group_drug), str(high_drug_mean), str(low_drug_mean), str(high_low_drug_diff), str(p_value)]) + "\n")
                    fout_log.write("Write data completed. " + str(datetime.now()) + "\n")

fout.close()
fout_log.close()
fout_error.close()

os.system("mv " + out_file + " " + out_folder)
os.system("mv " + out_log + " " + out_folder)
os.system("mv " + out_error + " " + out_folder)
print("move file completed.")
