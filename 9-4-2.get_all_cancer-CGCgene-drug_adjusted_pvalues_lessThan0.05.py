in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/CGC_lung19/mut_grouped_drug_prediction_t_test_p_values_CGC_all_cancer_adjusted_pvalues.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/CGC_lung19/mut_grouped_drug_prediction_t_test_p_values_CGC_all_cancer_adjusted_pvalues_lessThan0.05.txt"

import pandas as pd

df = pd.read_csv(in_file, sep = "\t")
df2 = df[df["adj_p_value"] < 0.05]
df2.to_csv(out_file, sep = "\t", index = False)
