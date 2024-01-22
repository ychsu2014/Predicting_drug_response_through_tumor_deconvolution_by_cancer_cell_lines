out_file = "4-9-2.predict_TCGA_cell_line_composition.sh"
model_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"
TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"

#scaden predict --model_dir /home/CBBI/hsuy1/projects/drugResponse/data/SCP542_UMI_processed/simulation_onlyNormalized/for_PRISM_by_cancertypes/processed_simulation_data/Breast_Cancer/model /home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/RNA/BC_RNA_geneSubset.txt --outname /home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/RNA/BC_RNA_geneSubset_prediction.txt

import os

simu_files = os.listdir(model_path)
model_folders = []
for simu_file in simu_files:
    temp = model_path + simu_file
    if os.path.isdir(temp) == True:
        model_folders.append(simu_file)

fout = open(out_file, "w")
for model_folder in model_folders:
    TCGA_file = model_folder.replace("_", "") + "_RNA_overlapped_genes_geneSubset.txt"
    TCGA_pred_file = TCGA_file.replace(".txt", "_prediction.txt")
    fout.write("scaden predict --model_dir " + model_path + model_folder + "/model " + TCGA_path + TCGA_file + " --outname " + TCGA_path + TCGA_pred_file + "\n")
fout.close()
