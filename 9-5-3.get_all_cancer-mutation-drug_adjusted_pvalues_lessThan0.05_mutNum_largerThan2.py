in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05.txt"
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05_mutNum_moreThan2_v3.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05_extremeSig_v3.txt"

up_diff_value = 0.5
down_diff_value = -0.5
sig_p_value = 10 ** (-50)

f = open(in_file)
fout1 = open(out_file1, "w")
fout2 = open(out_file2, "w")

lines = f.readlines()
header = lines[0]
fout1.write(header)
fout2.write(header)
lines = lines[1:]
count = 0
print(len(lines))
for line in lines:
    count += 1
    if count % 1000000 == 0:
        print(count)
    cols = line.strip("\n").split("\t")
    gene = cols[1]
    pro_mut = cols[2]
    drug = cols[3]
    mut_num = int(cols[4])
    mut_diff = float(cols[8])
    adj_p = float(cols[10])
    if gene != "" and pro_mut != "" and drug != "":
        if mut_num > 2:
            fout1.write(line)
        if mut_diff > up_diff_value and adj_p < sig_p_value:
            fout2.write(line)
        if mut_diff < down_diff_value and adj_p < sig_p_value:
            fout2.write(line)

f.close()
fout1.close()
fout2.close()
