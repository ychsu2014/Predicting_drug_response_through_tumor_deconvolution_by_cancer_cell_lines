#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_combNum.txt"
in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05_combNum.txt"


f = open(in_file)
cancer_totalNum_dict = {}
count = 0
for line in f:
    count += 1
    if count % 1000000 == 0:
        print(count)
    cols = line.strip("\n").split("\t")
    if cols[1] != "" and cols[2] != "" and cols[3] != "":
        cancer_type = cols[0]
        if cancer_type in cancer_totalNum_dict.keys():
            old_num = cancer_totalNum_dict[cancer_type]
            old_num += 1
            cancer_totalNum_dict[cancer_type] = old_num
        else:
            cancer_totalNum_dict[cancer_type] = 1

fout = open(out_file, "w")
for cancer_type in cancer_totalNum_dict.keys():
    total_num = cancer_totalNum_dict[cancer_type]
    fout.write(cancer_type + "\t" + str(total_num) + "\n")

f.close()
fout.close()
