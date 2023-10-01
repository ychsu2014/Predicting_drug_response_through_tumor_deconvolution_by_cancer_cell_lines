#dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized"
#simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data"
### simulation
#scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath

dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/lung_cancer_all_oncoKB_19"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/lung_cancer_19"
### simulation
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath
