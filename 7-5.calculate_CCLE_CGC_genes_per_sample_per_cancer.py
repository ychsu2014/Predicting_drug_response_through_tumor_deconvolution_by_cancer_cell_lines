in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered_194_CGC_matrix.txt"
in_cancer = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered_194_CGC_gene_count_per_cancer_type.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered_194_CGC_gene_count_per_sample.txt"

import pandas as pd
import numpy as np

# cancer type to list of sIDs
# sID to cancer type
f = open(in_cancer)
lines = f.readlines()
lines = lines[1:]
cancerType_sID_dict = {}
sID_cancerType_dict = {}
sID_cell_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    sID = cols[4]
    cancer_type = cols[1]
    cell = cols[3]
    sID_cancerType_dict[sID] = cancer_type
    sID_cell_dict[sID] = cell
    if cancer_type in cancerType_sID_dict.keys():
        old_list = cancerType_sID_dict[cancer_type]
        old_list.append(sID)
        old_list = list(set(old_list))
        cancerType_sID_dict[cancer_type] = old_list
    else:
        cancerType_sID_dict[cancer_type] = [sID]
f.close()

# cancer type to gene list
df = pd.read_csv(in_file, sep = "\t", index_col = 0)
cancerType_genes_dict = {}
for cancer_type in cancerType_sID_dict.keys():
    sID_list = cancerType_sID_dict[cancer_type]
    gene_list = df.columns[np.where(df.loc[sID_list].sum())]
    cancerType_genes_dict[cancer_type] = gene_list

# write CGC list per cancer
fout1 = open(out_file1, "w")
for cancer_type in cancerType_genes_dict.keys():
    gene_list = cancerType_genes_dict[cancer_type]
    gene_num = len(gene_list)
    fout1.write(cancer_type + "\t" + str(gene_num) + "\t" + "\t".join(gene_list) + "\n")
fout1.close()

# write CGC list per sample
fout2 = open(out_file2, "w")
sID_list = list(df.index)
for sID in sID_list:
    sID_gene_list = df.columns[np.where(df.loc[sID])]
    sID_gene_num = len(sID_gene_list)
    cancer_type = sID_cancerType_dict[sID]
    cell_name = sID_cell_dict[sID]
    fout2.write(cell_name + "\t" + cancer_type + "\t" + str(sID_gene_num) + "\t" + "\t".join(sID_gene_list) + "\n")
fout2.close()
