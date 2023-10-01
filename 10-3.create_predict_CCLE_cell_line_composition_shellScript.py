out_file = "10-4.predict_CCLE_cell_line_composition.sh"
model_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"
in_lung_model_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/processed_simulation_data/"
CCLE_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"

#scaden predict --model_dir /home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/for_PRISM_by_cancertypes/processed_simulation_data/Breast_Cancer/model /home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/RNA/BC_RNA_geneSubset.txt --outname /home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/RNA/BC_RNA_geneSubset_prediction.txt

import os

all_files = os.listdir(model_path)
model_folders = []
for afile in all_files:
    if "prediction" in afile:
        cancer_type = afile.replace("_prediction.txt", "")
        if cancer_type == "Lung_Cancer":
            temp = in_lung_model_path + cancer_type
        else:
            temp = model_path + cancer_type
        if os.path.isdir(temp) == True:
            model_folders.append(cancer_type)
        
fout = open(out_file, "w")
for model_folder in model_folders:
    CCLE_file = model_folder.replace("_", "") + "_RNA_ovelapped_genes_geneSubset.txt"
    CCLE_pred_file = CCLE_file.replace(".txt", "_prediction.txt")
    if model_folder == "Lung_Cancer":
        fout.write("scaden predict --model_dir " + in_lung_model_path + model_folder + "/model " + CCLE_path + CCLE_file + " --outname " + CCLE_path + CCLE_pred_file + "\n")
    else:
        fout.write("scaden predict --model_dir " + model_path + model_folder + "/model " + CCLE_path + CCLE_file + " --outname " + CCLE_path + CCLE_pred_file + "\n")
fout.close()
