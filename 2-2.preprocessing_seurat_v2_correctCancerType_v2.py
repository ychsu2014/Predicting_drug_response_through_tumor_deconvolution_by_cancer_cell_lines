# in Metadata.txt, only filtered cells were kept.
# in the UMI file, all the cells were kept.
# used 2-1.ID_cellLine_from_UMI_header.py to combine the cell line information in Metadata.txt and the UMI file.
# the output of 2-1 is "ID_cellline_cancertype.txt". This file is then manually checked, and saved as "cellline_cancertype_manuallyChecked.txt"
in_data = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_files/"
in_cell_line = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/ID_cellline_cancertype_v2_manuallyChecked_noHeader.txt"
results_file = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/Seurat_all_steps_result_record/pan_cancer.h5ad"
out_normalized_data = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/pan_cancer_norm_counts_all.txt"
out_normalized_ID_num = "/home/yuching/projects/drugResponse/data/scRNA/preprocessing_onlyNormalized/pan_cancer_data/ID_num.txt"

import numpy as np
import pandas as pd
import scanpy as sc
import random

# for cell filtering
min_genes_threshold = 200
# for gene filtering
min_cells_threshold = 3
# for cell filtering
n_genes_by_counts_threshold = 6000
pct_counts_mt_threshold = 13

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# cell line, cancer type information
f = open(in_cell_line)
lines = f.readlines()
ID_cellLine_dict = {}
ID_cancerType_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    ID = cols[0]
    cell_line = cols[1]
    cancer_type = cols[2]
    ID_cellLine_dict[ID] = cell_line
    ID_cancerType_dict[ID] = cancer_type

cell_line_df = pd.DataFrame.from_dict(ID_cellLine_dict, orient = "index", columns = ["cell_line"])
cancer_type_df = pd.DataFrame.from_dict(ID_cancerType_dict, orient = "index", columns = ["cancer_type"])

adata = sc.read_10x_mtx(in_data, var_names = "gene_symbols", cache =  True)

adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes = min_genes_threshold)
sc.pp.filter_genes(adata, min_cells = min_cells_threshold)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True, save = "_violinPlots.png")

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save = "_totalCounts_pctCountsMt_scatter.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save = "_totalCounts_nGenesByCounts.png")

adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts_threshold, :]
adata = adata[adata.obs.pct_counts_mt < pct_counts_mt_threshold, :]

sc.pp.normalize_total(adata, target_sum=1e4)

# save only normalized data for Scaden
df = pd.DataFrame(adata.X.todense())
df.columns = adata.var.index
df.to_csv(out_normalized_data, sep="\t")

# save ID to num table
fout_ID = open(out_normalized_ID_num, "w")
ID_list = list(adata.obs.index)
ID_num = -1
for ID in ID_list:
    ID_num += 1
    fout_ID.write(str(ID_num) + "\t" + ID + "\n")

fout_ID.close()
print("Saved files needed for Scaden.")