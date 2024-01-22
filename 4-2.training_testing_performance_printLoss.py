simu_path = "/home/yuching/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/"
CCC_all_results = "/home/yuching/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/CCC_all_results.txt"
TCGA_path = "/home/yuching/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/RNA/by_cancertypes/"

import os

processed_simu_path = simu_path + "processed_simulation_data"

# for rerun skin cancer
#simu_file_list = ["Skin_Cancer.h5ad"]

simu_file_list = []
for afile in os.listdir(simu_path):
    if ".h5ad" in afile:
        simu_file_list.append(afile)

# for testing
#simu_file_list = ["Colon_Colorectal_Cancer.h5ad"]

# two types of cancers were not in TCGA
simu_file_list.remove("Bone_Cancer.h5ad")
simu_file_list.remove("Neuroblastoma.h5ad")

for simufile in simu_file_list:
    print(simufile)
    TCGA_txt_file = TCGA_path + simufile.replace("_", "").replace(".h5ad", "_RNA_overlapped_genes.txt")
    # split to training/testing data
    os.system("python3 4-.split_to_training_testing.py " + simu_path + " " + processed_simu_path + " " + simufile)
    # reformat testing data from h5ad file to txt file
    #os.system("python3 4-4.reformat_bulk_RNA_data.py " + processed_simu_path + " " + simufile)
    # testing data processing
    test_h5ad_file = processed_simu_path + "/testing_" + simufile
    #test_txt_file = processed_simu_path + "/testing_" + simufile[:-5] + ".txt"
    test_processed_file = processed_simu_path + "/testing_" + simufile[:-5] + "_processed.h5ad"
    # scaden process $simuPath/testing_pan_cancer.h5ad $simuPath/testing_pan_cancer.txt --processed_path $simuPath/testing_pan_cancer_processed.h5ad
    #os.system("scaden process " + test_h5ad_file + " " + test_txt_file + " --processed_path " + test_processed_file)
    os.system("scaden process " + test_h5ad_file + " " + TCGA_txt_file + " --processed_path " + test_processed_file)
    # training data processing
    train_h5ad_file = processed_simu_path + "/training_" + simufile
    train_processed_file = processed_simu_path + "/training_" + simufile[:-5] + "_processed.h5ad"
    #os.system("scaden process " + train_h5ad_file + " " + test_txt_file + " --processed_path " + train_processed_file)
    os.system("scaden process " + train_h5ad_file + " " + TCGA_txt_file + " --processed_path " + train_processed_file)
    # reformat processed testing data from h5ad file to txt file
    os.system("python3 4-4.reformat_testing_data.py " + processed_simu_path + " " + simufile)
    # training
    #scaden train $simuPath/training_pan_cancer_processed.h5ad --model_dir $simuPath/model --steps 5000
    model_folder = processed_simu_path + "/" + simufile[:-5]
    os.system("mkdir " + model_folder)
    model_folder = model_folder + "/model"
    out_loss_record = model_folder + "/" + simufile[:-5] +"_loss_record.txt"
    #os.system("scaden train " + train_processed_file + " --model_dir " + model_folder + " --steps 5000")
    os.system("python3 4-5.training_scaden_model_printLoss.py " + train_processed_file + " " + model_folder + " " + out_loss_record)
    # testing
    #scaden predict --model_dir $simuPath/model $simuPath/testing_pan_cancer_processed.txt --outname $simuPath/pan_cancer_prediction.txt
    test_processed_txt_file = processed_simu_path + "/testing_" + simufile[:-5] + "_processed.txt"
    predict_result = processed_simu_path + "/" + simufile[:-5] + "_prediction.txt"
    os.system("scaden predict --model_dir " + model_folder + " " + test_processed_txt_file + " --outname " + predict_result)
    # calculate concordance correlation coefficient
    os.system("python3 4-6.CCC_evaluation_pooledAllForOneDataset.py " + processed_simu_path + " " + simufile + " " + CCC_all_results)
    print("Completed!")
