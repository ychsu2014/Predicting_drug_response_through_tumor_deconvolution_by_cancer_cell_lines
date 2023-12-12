in_scRNA_lung = "/home/yuching/drugResponse/data/scRNA/preprocessing_onlyNormalized/by_cancertypes_all_cell_lines/Lung_Cancer_norm_counts_all.txt"
in_scRNA_lung_celltypes = "/home/yuching/drugResponse/data/scRNA/preprocessing_onlyNormalized/by_cancertypes_all_cell_lines/Lung_Cancer_celltypes.txt"
in_gene = "/home/yuching/drugResponse/data/scRNA/preprocessing_onlyNormalized/keep_only_TCGA_genes_lung_19/Lung_Cancer_norm_counts_all.txt"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import matplotlib

matplotlib.use('agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

sys.setrecursionlimit(100000000)

df = pd.read_csv(in_scRNA_lung, sep = "\t", index_col = 0)
cell_df = pd.read_csv(in_scRNA_lung_celltypes, sep = "\t", index_col = 0)
gene_df = pd.read_csv(in_gene, sep = "\t", index_col = 0)
df =  pd.concat([df, cell_df], axis = 1)
df.set_index("Celltype", inplace = True)

df = df[gene_df.columns]
print("dataframe completed.")

fig = sns.clustermap(df, vmax = 10)
fig.savefig("lung_heatmap_10.pdf", format = "pdf")
plt.close()
print("fig completed.")