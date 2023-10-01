dataPath="/home/yuching/projects/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/by_cancertypes_selected_cell_lines"
simuPath="/home/yuching/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized"

### simulation
scaden simulate -d $dataPath -c 500 -n 8000 --pattern "*_norm_counts_all.txt" -o $simuPath