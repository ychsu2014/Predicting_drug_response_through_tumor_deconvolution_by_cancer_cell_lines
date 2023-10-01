in_h5ad = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/Seurat_all_steps_result_record/pan_cancer.h5ad"
# the sample ID must be able to be converted to float.
out_celltype = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/pan_cancer_celltypes.txt"

import scanpy as sc

fout_celltype = open(out_celltype, "w")

pan_h5ad = sc.read_h5ad(in_h5ad)
gene_list = list(pan_h5ad.var["gene_ids"].index)
ID_list = list(pan_h5ad.obs.index)
cell_line_list = list(pan_h5ad.obs["cell_line"])

# write
fout_celltype.write("\tCelltype\n")
print(len(ID_list))
for i in range(len(ID_list)):
    ID_num = str(i)
    if i % 5000 == 0:
        print(i)
    fout_celltype.write(ID_num + "\t")
    fout_celltype.write(cell_line_list[i] + "\n")

fout_celltype.close()
