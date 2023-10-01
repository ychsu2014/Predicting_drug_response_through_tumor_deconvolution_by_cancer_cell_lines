in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues_lessThan_0.05.txt"

f = open(in_file)
fout = open(out_file, "w")

adj_pvalue_colnum = 10

lines = f.readlines()
header = lines[0]
fout.write(header)
lines = lines[1:]
count = 0
for line in lines:
    count += 1
    if count % 1000000 == 0:
        print(count)
    cols = line.strip("\n").split("\t")
    adj_pvalue = float(cols[adj_pvalue_colnum])
    if adj_pvalue < 0.05:
        fout.write(line)

f.close()
fout.close()
