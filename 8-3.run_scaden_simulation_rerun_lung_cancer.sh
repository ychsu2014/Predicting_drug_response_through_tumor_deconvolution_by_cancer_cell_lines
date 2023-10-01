mkdir /home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_oncoKB_with_drugs_11/
mkdir /home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/
mkdir /home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_mut_and_CGC_genes_30/

dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/lung_cancer_oncoKB_with_drugs_11/"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_oncoKB_with_drugs_11/"

### simulation
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath


dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/lung_cancer_all_oncoKB_19/"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/"

### simulation
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath


dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/lung_cancer_all_oncoKB_mut_and_CGC_genes_30/"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_mut_and_CGC_genes_30/"

### simulation
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath
