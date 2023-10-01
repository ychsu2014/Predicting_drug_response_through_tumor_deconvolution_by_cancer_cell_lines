in_file = "/home/yuching/projects/drugResponse/data/CCLE/RNA/CCLE_all_RNA.h5ad"
#idx stripped_cell_line_name cell_line_name  CCLE_Name
#0   LC1SQSF LC-1/sq-SF  LC1SQSF_LUNG
in_ID = "/home/yuching/projects/drugResponse/data/CCLE/RNA/CCLE_all_RNA_idx_depmapID.txt"
#scRNA   scRNA_cancerType    CCLE_full   CCLE_short  depmap_ID
#SW579_THYROID   Thyroid Cancer  SW579   SW579   ACH-000163
in_phe = "/home/yuching/projects/drugResponse/data/processed_drug_datasets/PRISM_scRNA_CCLE_cell_lines_list_v2_lung19.txt"
out_path1 = "/home/yuching/projects/drugResponse/data/CCLE/RNA/by_cancer/h5ad/"
out_path2 = "/home/yuching/projects/drugResponse/data/CCLE/RNA/by_cancer/idx_files/"

out1_suffix = "_RNA.h5ad"
out2_suffix = "_RNA_idx_sID.txt"

import anndata as ad

### the selected cell lines in the subsequent analysis
# scRNAcancerType to depmap ID
scRNACancerType_cellID_dict = {}
f = open(in_phe)
lines = f.readlines()
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    cancer_type = cols[1]
    depmap_ID = cols[4]
    if cancer_type in scRNACancerType_cellID_dict.keys():
        ID_list = scRNACancerType_cellID_dict[cancer_type]
        ID_list.append(depmap_ID)
        ID_list = list(set(ID_list))
        scRNACancerType_cellID_dict[cancer_type] = ID_list
    else:
        scRNACancerType_cellID_dict[cancer_type] = [depmap_ID]
f.close()

# depmap ID to idx lines
# only keep the depmap IDs included in the subsequent analysis
f = open(in_ID)
lines = f.readlines()
header = lines[0]
lines = lines[1:]
ID_idxLine_dict = {}
ID_idx_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    depmap_ID = cols[4]
    idx = cols[0]
    ID_idxLine_dict[depmap_ID] = line
    ID_idx_dict[depmap_ID] = idx
f.close()

# write idx ID mapping table by cancer
# key: cancer_type, value: list of idx
cancerType_idx_dict = {}
for cancer_type in scRNACancerType_cellID_dict.keys():
    depmap_ID_list = scRNACancerType_cellID_dict[cancer_type]
    cancer_type = cancer_type.replace(" ", "").replace("/", "")
    fout2 = open(out_path2 + cancer_type + out2_suffix, "w")
    idx_list = []
    for depmap_ID in depmap_ID_list:
        idx_line = ID_idxLine_dict[depmap_ID]
        fout2.write(idx_line)
        idx = ID_idx_dict[depmap_ID]
        idx_list.append(idx)
    cancerType_idx_dict[cancer_type] = idx_list 
    fout2.close()

all_RNA =  ad.read_h5ad(in_file)
for cancer_type in cancerType_idx_dict.keys():
    print(cancer_type)
    idx_list = cancerType_idx_dict[cancer_type]
    RNA_temp = all_RNA[idx_list]
    RNA_temp.write_h5ad(out_path1 + cancer_type + out1_suffix)