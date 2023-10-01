dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/preprocessing_onlyNormalized/oncoKB_CGC"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_CGC"
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath

dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/preprocessing_onlyNormalized/oncoKB"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB"
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath

dataPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/preprocessing_onlyNormalized/oncoKB_actionable_target"
simuPath="/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_actionable_target"
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath
