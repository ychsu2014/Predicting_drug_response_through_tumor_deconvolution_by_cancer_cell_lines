#in_TCGA = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7861_CGC_matrix.txt"
#in_CCLE = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered_194_CGC_matrix.txt"
#in_CCLE_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
#in_TCGA_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7861.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_CGC/TCGA_CCLE_overlapped_CGC_genes.txt"

in_TCGA = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7781_CGC_matrix.txt"
in_CCLE = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194_filteredLungCancer19_CGC_matrix.txt"
in_CCLE_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list.txt"
in_TCGA_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7781.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_CGC/TCGA_CCLE_overlapped_CGC_genes.txt"
#out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_CGC/TCGA_CCLE_overlapped_CGC_genes_per_sample.txt"
in_lung_ID = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/selected_cell_lines_for_lung_cancer_all_oncoKB_19.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_CGC/lung_cancer_all_oncoKB_19/TCGA_CCLE_overlapped_CGC_genes_all_oncoKB_19.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_CGC/lung_cancer_all_oncoKB_19/TCGA_CCLE_overlapped_CGC_genes_per_sample_all_oncoKB_19.txt"

import pandas as pd
import numpy as np

TCGA = pd.read_csv(in_TCGA, sep = "\t", index_col = 0)
CCLE = pd.read_csv(in_CCLE, sep = "\t", index_col = 0)

# rerun lung cancer
f = open(in_lung_ID)
lines = f.readlines()
lung_19_ID = []
for line in lines:
    lung_19_ID.append(line.strip("\n"))
f.close()

def parseToDict(inFile, inCancerColnum, inIDColnum):
    f = open(inFile)
    lines = f.readlines()
    lines = lines[1:]
    outCancerIDDict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        cancerType = cols[inCancerColnum]
        sID = cols[inIDColnum]
        if cancerType in outCancerIDDict.keys():
            oldList = outCancerIDDict[cancerType]
            oldList.append(sID)
            oldList = list(set(oldList))
            outCancerIDDict[cancerType] = oldList
        else:
            outCancerIDDict[cancerType] = [sID]
    f.close()
    return(outCancerIDDict)

# key: scRNA cancer types, value: depmap IDs
CCLE_cancerType_ID_dict = parseToDict(in_CCLE_ID, 1, 4)

# for rerun lung cancer
CCLE_cancerType_ID_dict["Lung Cancer"] = lung_19_ID

# key: scRNA cancer types, value: TCGA sIDs
TCGA_cancerType_ID_dict = parseToDict(in_TCGA_ID, 2, 0)

fout = open(out_file, "w")
fout2 = open(out_file2, "w")
fout.write("cancer_type\tgene_num\tgenes\n")
fout2.write("TCGA_ID\tcancer_type\tgene_num\tgenes\n")
for cancer_type in TCGA_cancerType_ID_dict:
    TCGA_ID_list = TCGA_cancerType_ID_dict[cancer_type]
    TCGA_temp_idx = np.where(TCGA.loc[TCGA_ID_list].sum())
    TCGA_genes = TCGA.columns[TCGA_temp_idx]
    CCLE_ID_list = CCLE_cancerType_ID_dict[cancer_type]
    CCLE_temp_idx = np.where(CCLE.loc[CCLE_ID_list].sum())
    CCLE_genes = CCLE.columns[CCLE_temp_idx]
    TCGA_CCLE_overlapped_genes = list(set(TCGA_genes).intersection(set(CCLE_genes)))
    fout.write(cancer_type)
    fout.write("\t" + str(len(TCGA_CCLE_overlapped_genes)))
    for overlapped_gene in TCGA_CCLE_overlapped_genes:
        fout.write("\t" + overlapped_gene)
    fout.write("\n")
    # write to individual TCGA samples
    for TCGA_ID in TCGA_ID_list:
        TCGA_s_temp_idx = np.where(TCGA.loc[TCGA_ID] != 0)
        TCGA_s_genes = TCGA.columns[TCGA_s_temp_idx]
        TCGA_s_CCLE_genes = list(set(TCGA_s_genes).intersection(set(CCLE_genes)))
        fout2.write(TCGA_ID + "\t" + cancer_type + "\t" + str(len(TCGA_s_CCLE_genes)))
        for agene in TCGA_s_CCLE_genes:
            fout2.write("\t" + agene)
        fout2.write("\n")

fout.close()
fout2.close()
