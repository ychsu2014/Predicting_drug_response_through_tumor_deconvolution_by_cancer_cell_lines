in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues.txt"
file_prefix = "mut_grouped_drug_prediction_t_test_p_values_all_mutations_"
pvalue_colnum = 9

import os
from statsmodels.stats.multitest import multipletests

# t-test results for all CGC genes
all_files = os.listdir(in_path)
all_mutation_file_list = []
for afile in all_files:
    if file_prefix in afile and "all_cancer" not in afile:
        all_mutation_file_list.append(afile)        

# combine all the data
all_data_lines = []
for amutfile in all_mutation_file_list:
    f = open(in_path + amutfile)
    lines = f.readlines()
    lines = lines[1:]
    all_data_lines += lines
    f.close()
print("All data combined.")

# check num of cols are the same
column_num_check = len(all_data_lines[0].strip("\n").split("\t"))
for adataline in all_data_lines:
    cols = adataline.strip("\n").split("\t")
    if len(cols) != column_num_check:
        print("Please check!!!!!!")
        print(adataline)

# all p-values list
p_value_list = []
for data_line in all_data_lines:
    cols = data_line.strip("\n").split("\t")
    p_value = float(cols[pvalue_colnum])
    p_value_list.append(p_value)
print("Number of p-values: " + str(len(p_value_list)))

# bonferroni correction
adj_pvalues = list(multipletests(p_value_list, method = "bonferroni")[1])
print("Number of all data lines: " + str(len(all_data_lines)))
print("Number of adjusted p-values: " + str(len(adj_pvalues)))
print("Bonferroni correction completed.")

# write to file
fout = open(out_file, "w")
fout.write("\t".join(["Cancer_type", "gene", "protein_hgvs", "drug", "mut_num", "wt_num", "mut_mean", "wt_mean", "mut_wt_diff", "p_values", "adj_p_value"]) + "\n")
count = 0
for i in range(len(all_data_lines)):
    count += 1
    if count % 100000 == 0:
        print(count)
    data_line = all_data_lines[i]
    adj_pvalue = adj_pvalues[i]
    fout.write(data_line.strip("\n") + "\t" + str(adj_pvalue) + "\n")
fout.close()
