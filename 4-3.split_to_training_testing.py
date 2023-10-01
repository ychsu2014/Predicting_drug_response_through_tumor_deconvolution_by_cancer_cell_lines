#in_simu = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/pan_cancer.h5ad"
#out_train = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/training_pan_cancer.h5ad"
#out_test = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/testing_pan_cancer.h5ad"

#in_simu = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/by_cancerType/Bladder_Cancer.h5ad"
#out_train = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/by_cancerType/training_Bladder_Cancer.h5ad"
#out_test = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation/by_cancerType/testing_Bladder_Cancer.h5ad"

#in_simu = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/Breast_Cancer.h5ad"
#out_train = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/training_Breast_Cancer.h5ad"
#out_test = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/testing_Breast_Cancer.h5ad"

#in_simu = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/Colon_Colorectal_Cancer.h5ad"
#out_train = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/training_Colon_Colorectal_Cancer.h5ad"
#out_test = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/testing_Colon_Colorectal_Cancer.h5ad"

#in_simu = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/Head_and_Neck_Cancer.h5ad"
#out_train = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/training_Head_and_Neck_Cancer.h5ad"
#out_test = "/home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/by_cancertypes/testing_Head_and_Neck_Cancer.h5ad"

# input file path to simulation data & the filename of the simulation data 

import scanpy as sc
import random
import sys


data_path = sys.argv[1]
out_path = sys.argv[2]
file_suffix = sys.argv[3]

in_simu = data_path + "/" + file_suffix
out_train = out_path + "/training_" + file_suffix
out_test = out_path + "/testing_" + file_suffix

# indexes for testing data
random.seed(0)
test_idx = random.sample(range(0,7999), 1600)
print(len(test_idx))

# indexes for training data
all_idx = list(range(8000))
train_idx = list(set(all_idx) - set(test_idx))
print(len(train_idx))

pan_cancer = sc.read_h5ad(in_simu)
train_pan = pan_cancer[pan_cancer.obs_names[train_idx]]
test_pan = pan_cancer[pan_cancer.obs_names[test_idx]]
train_pan.write(out_train)
test_pan.write(out_test)
