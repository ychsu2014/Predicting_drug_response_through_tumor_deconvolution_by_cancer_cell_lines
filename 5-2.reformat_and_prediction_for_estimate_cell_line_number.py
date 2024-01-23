#simu_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/"
#out_result = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/CCC_all_results.txt"
TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/h5ad/"
model_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"

simu_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/lung_cancer_19/"
out_result = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/lung_cancer_19/CCC_lung_cancer_19_results.txt"
lung_model_folder = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/processed_simulation_data/Lung_Cancer/model/"

import os

all_files = os.listdir(simu_path)
simu_file_list = []
for afile in all_files:
    # two types of cancers were not in TCGA
    if "h5ad" in afile and "processed" not in afile and "Bone_Cancer" not in afile and "Neuroblastoma" not in afile:
        simu_file_list.append(afile)

for simufile in simu_file_list:
    print(simufile)
    TCGA_txt_file = TCGA_path + simufile.replace("_", "").replace(".h5ad", "_RNA.txt")
    if simufile[:-5] == "Lung_Cancer":
        model_folder = lung_model_folder
    else:
        model_folder = model_path + simufile[:-5] + "/model"
    # data processing
    simufile_full = simu_path + simufile
    processed_h5ad = simufile[:-5] + "_processed.h5ad"
    processed_h5ad_full = simu_path + simufile[:-5] + "_processed.h5ad"
    os.system("scaden process " + simufile_full + " " + TCGA_txt_file + " --processed_path " + processed_h5ad_full)
    # data reformat from h5ad to txt
    os.system("python3 11-2_1.reformat_h5ad_to_txt.py " + simu_path + " " + processed_h5ad)
    processed_txt = simu_path + simufile[:-5] + "_processed.txt"
    predict_result = simu_path + simufile[:-5] + "_prediction.txt"
    # prediction
    os.system("scaden predict --model_dir " + model_folder + " " + processed_txt + " --outname " + predict_result)
    # calculate concordance correlation coefficient
    os.system("python3 11-2_2.CCC_evaluation_pooledAllForOneDataset.py " + simu_path + " " + simufile + " " + out_result)
