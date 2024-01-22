out_file = "4-9-2-1.predict_CCLE_cell_line_composition.sh"
model_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"
CCLE_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"

import os

simu_files = os.listdir(model_path)
model_folders = []
for simu_file in simu_files:
    temp = model_path + simu_file
    if os.path.isdir(temp) == True:
        model_folders.append(simu_file)

fout = open(out_file, "w")
for model_folder in model_folders:
    CCLE_file = model_folder.replace("_", "") + "_RNA_overlapped_genes_geneSubset.txt"
    CCLE_pred_file = CCLE_file.replace(".txt", "_prediction.txt")
    fout.write("scaden predict --model_dir " + model_path + model_folder + "/model " + CCLE_path + CCLE_file + " --outname " + CCLE_path + CCLE_pred_file + "\n")
fout.close()
