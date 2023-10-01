# the output data needs to be normalized and log2-transformed.
# need to test this script by the previous simulation dataset (done)

# prediction files of all cell lines
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"
in_lung_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/lung_cancer_all_oncoKB_19/"
# scRNA count data
in_scRNA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/keep_only_TCGA_genes_lung_19/"
# pseudo bulk RNA data
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/pseudo_bulk_RNA/"

prediction_suffix = "_RNA_overlapped_genes_geneSubset_prediction.txt"

import numpy as np
import pandas as pd
import os
import anndata as ad

# scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath
sample_size = 500
random_seed_num = 100
#num_samples = 8000

all_files = os.listdir(in_scRNA_path)
count_files = []
for afile in all_files:
    if "_norm_counts_all.txt" in afile:
        count_files.append(afile)

# for debugging
count_files = ["Lung_Cancer_norm_counts_all.txt"]

# bone cancer, neuroblastoma are not in TCGA cancer types
if "Neuroblastoma_norm_counts_all.txt" in count_files:
    count_files.remove("Neuroblastoma_norm_counts_all.txt")
if "Bone_Cancer_norm_counts_all.txt" in count_files:
    count_files.remove("Bone_Cancer_norm_counts_all.txt")


#Options:
#  -o, --out TEXT           Directory to store output files in
#  -d, --data TEXT          Path to scRNA-seq dataset(s)
#  -c, --cells INTEGER      Number of cells per sample [default: 100]
#  -n, --n_samples INTEGER  Number of samples to simulate [default: 1000]
#  --pattern TEXT           File pattern to recognize your processed scRNA-seq
#                           count files
#  -u, --unknown TEXT       Specifiy cell types to merge into the unknown
#                           category. Specify this flag for every cell type you
#                           want to merge in unknown. [default: unknown]
#  -p, --prefix TEXT        Prefix to append to training .h5ad file [default:
#                           data]
#  -f, --data-format TEXT   Data format of scRNA-seq data, can be 'txt' or
#                           'h5ad' [default: 'txt']
#  --help                   Show this message and exit.

for count_file in count_files:
    cancer_type = count_file.replace("_norm_counts_all.txt", "")
    print(cancer_type)
    cell_type_file = count_file.replace("_norm_counts_all.txt", "_celltypes.txt")
    pred_file = count_file.replace("_norm_counts_all.txt", "").replace("_", "") + prediction_suffix
    # load scRNA dataset
    data_x = pd.read_table(in_scRNA_path + count_file, index_col = 0)
    data_y = pd.read_table(in_scRNA_path + cell_type_file)
    celltypes = list(set(data_y["Celltype"].tolist()))
    # create subsample dataset
    no_avail_cts = len(celltypes)
    # use lung cancer 19 cell lines
    if cancer_type != "Lung_Cancer":
        f = open(in_path + pred_file)
    else:
        f = open(in_lung_path + pred_file)
    lines = f.readlines()
    header = lines[0].strip("\n").split("\t")[1:]
    lines = lines[1:]
    sim_x = []
    sim_y = []
    count = -1
    # key: new index, value: cell line name
    idx_sID_dict = {}
    for line in lines:
        count += 1
        cols = line.strip("\n").split("\t")
        idx_sID_dict[count] = cols[0]
        ## fracs: proportions for each cell line
        fracs = cols[1:]
        fracs = list(map(float, fracs))
        ## samp_fracs: number of cells to be selected for each cell line
        samp_fracs = np.multiply(fracs, sample_size)
        samp_fracs = list(map(int, samp_fracs))
        #print(cols[0])
        #print(samp_fracs)
        # available_celltypes = celltypes
        # fracs_complete = fracs
        artificial_samples = []
        for i in range(no_avail_cts):
            ct = celltypes[i]
            ## cells_sub: cells x genes
            cells_sub = data_x.loc[np.array(data_y["Celltype"] == ct), :]    
            np.random.seed(seed = random_seed_num)
            cells_fraction = np.random.randint(0, cells_sub.shape[0], samp_fracs[i])
            cells_sub = cells_sub.iloc[cells_fraction, :]
            ## list of dataframes of the selected cells
            artificial_samples.append(cells_sub)
        ## dataframe of all the selected cells
        df_samp = pd.concat(artificial_samples, axis=0)
        ## pseudo bulk data
        df_samp = df_samp.sum(axis=0)
        sim_x.append(df_samp)
        sim_y.append(fracs)
    sim_x = pd.concat(sim_x, axis = 1).T
    ## sim_y, ratios: proportions of cell types
    sim_y = pd.DataFrame(sim_y, columns = header)
    # tmp_x = sim_x, tmp_y = sim_y
    sim_x = sim_x.sort_index(axis = 1)
    ratios = pd.DataFrame(sim_y, columns = header)
    ratios["ds"] = pd.Series(np.repeat(cancer_type, sim_y.shape[0]), index = ratios.index)
    ratios["sample_ID"] = pd.Series(idx_sID_dict)
    X = sim_x.to_numpy()
    ann_data = ad.AnnData(
        X,
        obs = ratios,
        var = pd.DataFrame(columns = [], index = list(sim_x)),
        dtype = X.dtype
    )
    ann_data.write(out_path + cancer_type + "_pseudo_bulk.h5ad")
