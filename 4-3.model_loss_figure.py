### for all cell line
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/processed_simulation_data/loss_figures/"

### for lung cancer 19
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/processed_simulation_data/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/scRNA/simulation_onlyNormalized/lung_cancer_all_oncoKB_19/processed_simulation_data/loss_figures/"

##### run after 12.* scripts #####
### for CGC_ oncoKB
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_CGC/processed_simulation_data/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_CGC/processed_simulation_data/loss_figures/"
### for all oncoKB except for lung cancer
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB/processed_simulation_data/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB/processed_simulation_data/loss_figures/"
### for oncoKB actionable targets
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_actionable_target/processed_simulation_data/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/greedy_selected_cell_line/simulation_onlyNormalized/oncoKB_actionable_target/processed_simulation_data/loss_figures/"
##########

import os
import matplotlib.pyplot as plt

all_files = os.listdir(in_path)
pred_files = []
for afile in all_files:
    if "_prediction.txt" in afile:
        pred_files.append(afile)


for pred_file in pred_files:
    cancer_type = pred_file.replace("_prediction.txt", "")
    model_path = in_path + cancer_type + "/model/"
    loss_file = model_path + cancer_type + "_loss_record.txt"
    f = open(loss_file)
    lines = f.readlines()
    m256_dict = {}
    m512_dict = {}
    m1024_dict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        model = cols[0]
        step = int(cols[1])
        loss = float(cols[3])
        if model == "256":
            m256_dict[step] = loss
        elif model == "512":
            m512_dict[step] = loss
        elif model == "1024":
            m1024_dict[step] = loss
    x_list = []
    m256_y_list = []
    m512_y_list = []
    m1024_y_list = []
    for i in range(1, len(m256_dict.keys()) +1):
        x_list.append(i)
        m256_loss = m256_dict[i]
        m256_y_list.append(m256_loss)
        m512_loss = m512_dict[i]
        m512_y_list.append(m512_loss)
        m1024_loss = m1024_dict[i]
        m1024_y_list.append(m1024_loss)
    line1, = plt.plot(x_list, m256_y_list, "b")
    line2, = plt.plot(x_list, m512_y_list, "r")
    line3, = plt.plot(x_list, m1024_y_list, "g")
    plt.legend(handles = [line1, line2, line3], labels = ["m256", "m512", "m1024"], loc = "best")
    plt.xlabel("Step")
    plt.ylabel("Loss")
    plt.title(cancer_type)
    plt.savefig(out_path + cancer_type + "_loss_plot.pdf", format = "pdf")
    plt.close()
