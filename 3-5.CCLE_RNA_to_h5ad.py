### RNA data of CCLE (depmap 2022Q2): log2(TPM+1)
in_RNA = "/home/yuching/projects/drugResponse/data/CCLE/RNA/CCLE_expression.csv"
in_sample = "/home/yuching/projects/drugResponse/data/CCLE/RNA/sample_info.csv"
out_h5ad = "/home/CBBI/hsuy1/projects/drugResponse/data/CCLE/RNA/CCLE_all_RNA.h5ad"
### idx, depmapID, stripped_cell_line_name (ex. SLR21), cell_line_name (ex.SLR 21), CCLE_Name (ex. SLR21_KIDNEY)
out_idx = "/home/CBBI/hsuy1/projects/drugResponse/data/CCLE/RNA/CCLE_all_RNA_idx_depmapID.txt"

from scipy.sparse import csr_matrix
import numpy as np
import anndata as ad

### depmapID to cell line names
f = open(in_sample)
lines = f.readlines()
lines = lines[1:]
depmapID_cellnames_dict = {}
for line in lines:
    cols = line.strip("\n").split(",")
    ID = cols[0]
    full_name = cols[1]
    stripped_name = cols[2]
    stripped_name_cancer = cols[3]
    depmapID_cellnames_dict[ID] = (stripped_name, full_name, stripped_name_cancer)    
f.close()

f = open(in_RNA)
lines = f.readlines()
header = lines[0]
### gene list
# TSPAN6 (7105)
gene_list = header.strip("\n").split(",")[1:]
gene_name_list = []
for gene in gene_list:
    gene_name = gene.split(" ")[0]
    gene_name_list.append(gene_name)
### write idx to cell line names mapping table
fout = open(out_idx, "w")
fout.write("idx\tstripped_cell_line_name\tcell_line_name\tCCLE_Name\tdepmap_ID\n")
lines = lines[1:]
count = -1
for line in lines:
    count += 1
    cols = line.strip("\n").split(",")
    ID = cols[0]
    if ID in depmapID_cellnames_dict.keys():
        (stripped_name, full_name, stripped_name_cancer) = depmapID_cellnames_dict[ID]
    else:
        stripped_name = ""
        full_name = ""
        stripped_name_cancer = ""
    fout.write("\t".join([str(count), stripped_name, full_name, stripped_name_cancer, ID]) + "\n")
fout.close()
f.close()

### prepare data for converting to compressed sparse row matrix (csr matrix)
# column: genes, row: cell lines
# pay attention to special value handling!!!
col_idx_list = []
row_idx_list = []
data_list = []
row_idx = -1
for line in lines:
    row_idx += 1
    cols = line.strip("\n").split(",")
    cols = cols[1:]
    col_idx = -1
    for col in cols:
        col_idx += 1
        col_idx_list.append(col_idx)
        row_idx_list.append(row_idx)
        data_list.append(float(col))

row_idx_array = np.array(row_idx_list)
col_idx_array = np.array(col_idx_list)
adata = csr_matrix((data_list, (row_idx_array, col_idx_array)), shape = (max(row_idx_array)+1, max(col_idx_array)+1))
adata = ad.AnnData(adata, dtype = float)
adata.var_names = gene_name_list

# write all RNA data to file
adata.write_h5ad(out_h5ad)
