#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/testing_pan_cancer_processed.h5ad"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/testing_pan_cancer_processed.txt"

#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/by_cancerType/test_run_bladder/testing_Bladder_Cancer_processed.h5ad"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/by_cancerType/test_run_bladder/testing_Bladder_Cancer_processed.txt"

#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/testing_Breast_Cancer_processed.h5ad"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/testing_Breast_Cancer_processed.txt"

#data_path = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/"
#in_file = data_path + "testing_Head_and_Neck_Cancer.h5ad"
#out_file = data_path + "testing_Head_and_Neck_Cancer.txt"

# input the file path to simulation data & the filename of the simulation data

import scanpy as sc
import pandas as pd
import sys

data_path = sys.argv[1]
simu_name = sys.argv[2]
in_file = data_path + "/testing_" + simu_name[:-5] + "_processed.h5ad"
out_file = data_path + "/testing_" + simu_name[:-5] + "_processed.txt"

# for log-transformed and scaled testing data
adata = sc.read_h5ad(in_file)
df = adata.to_df().T
df.to_csv(out_file, sep = "\t")
