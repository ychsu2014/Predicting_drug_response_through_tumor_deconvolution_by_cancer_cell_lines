in_path = "C:/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/"
in_idx_path = "C:/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/idx_files/"

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# bar plot & pie chart

all_files = os.listdir(in_path)
prediction_files = []
for afile in all_files:
    if "_RNA_geneSubset_addZero_prediction.txt" in afile:
        prediction_files.append(afile)

cancer_type_dict = {}
cancer_type_dict["BileDuctCancer"] = "Bile Duct Cancer"
cancer_type_dict["BladderCancer"] = "Bladder Cancer"
cancer_type_dict["BrainCancer"] = "Brain Cancer"
cancer_type_dict["BreastCancer"] = "Breast Cancer"
cancer_type_dict["ColonColorectalCancer"] = "Colon/Colorectal Cancer"
cancer_type_dict["EndometrialUterineCancer"] = "Endometrial/Uterine Cancer"
cancer_type_dict["EsophagealCancer"] = "Esophageal Cancer"
cancer_type_dict["GastricCancer"] = "Gastric Cancer"
cancer_type_dict["HeadandNeckCancer"] = "Head and Neck Cancer"
cancer_type_dict["KidneyCancer"] = "Kidney Cancer"
cancer_type_dict["LiverCancer"] = "Liver Cancer"
cancer_type_dict["LungCancer"] = "Lung Cancer"
cancer_type_dict["OvarianCancer"] = "Ovarian Cancer"
cancer_type_dict["PancreaticCancer"] = "Pancreatic Cancer"
cancer_type_dict["ProstateCancer"] = "Prostate Cancer"
cancer_type_dict["Sarcoma"] = "Sarcoma"
cancer_type_dict["SkinCancer"] = "Skin Cancer"
cancer_type_dict["ThyroidCancer"] = "Thyroid Cancer"

color_list = ["red", "orange", "limegreen", "turquoise", "gold", "blueviolet", "peru", "sienna", "olive", "magenta", "palevioletred", "tan", "lightgreen", "deeppink", "grey", "darkviolet", "darkgreen", "royalblue", "navy", "rosybrown"]

# key: cell, value: self proportion
cell_selfValue_dict = {}
# key: cancer, value: [cell names]
cancer_cell_dict = {}
for predict_file in prediction_files:
    cancer_title = cancer_type_dict[predict_file.replace("_RNA_geneSubset_addZero_prediction.txt", "")]
    idx_file = predict_file.replace("_RNA_geneSubset_addZero_prediction.txt", "_RNA_idx_sID.txt")
    idx_file = in_idx_path + idx_file
    # key: idx, value: cell name
    f = open(idx_file)
    lines = f.readlines()
    idx_dict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        idx_dict[cols[0]] = cols[1]
    predict_file = in_path + predict_file
    pred_df = pd.read_csv(predict_file, sep = "\t", index_col = 0)
    pred_df_idx_list = list(pred_df.index)
    # pie charts
    cell_name_list = []
    for aidx in pred_df_idx_list:
        plt.figure(figsize=(10, 10))
        mylabels = list(pred_df.loc[aidx].index)
        value_list = np.array(list(pred_df.loc[aidx]))
        plt.pie(value_list, labels = mylabels, textprops = {"fontsize":14})
        cell_name = idx_dict[str(aidx)]
        out_png = predict_file.replace("_RNA_geneSubset_addZero_prediction.txt","_" + cell_name + "_pie_chart.png")
        plt.title(cell_name, fontsize = 14)
        plt.tight_layout()
        plt.savefig(out_png)
        plt.close()
        # for stacked bar plot
        cell_name_list.append(cell_name)
        # for all cancer bar plot
        self_value = pred_df.loc[aidx][cell_name]
        cell_selfValue_dict[cell_name] = self_value
        if cancer_title in cancer_cell_dict.keys():
            old_list = cancer_cell_dict[cancer_title]
            old_list.append(cell_name)
            old_list = list(set(old_list))
            cancer_cell_dict[cancer_title] = old_list
        else:
            cancer_cell_dict[cancer_title] = [cell_name]
    # stacked bar plot per cancer
    # key: cell names, value: proportions
    stacked_bars = {}
    for col in pred_df.columns:
        stacked_bars[col] = np.array(pred_df[col])
    width = 0.5
    fig, ax = plt.subplots(figsize = (10,7))
    row_num = len(cell_name_list)
    bottom = np.zeros(row_num)
    # color list
    color_count = -1
    for category, count in stacked_bars.items():
        color_count += 1
        bar_color = color_list[color_count]
        p = ax.bar(cell_name_list, count, width, label = category, bottom = bottom, color = bar_color)
        bottom += count
    ax.set_title(cancer_title)
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    plt.xticks(rotation = 45)
    #ax.legend(loc="upper right")
    out_bar = predict_file.replace("_RNA_geneSubset_addZero_prediction.txt","") + "_stacked_bar_plot.png"
    plt.tight_layout()
    plt.savefig(out_bar)
    plt.close()

# prepare the data for all cancer stacked bar plot
self_pro_list = []
other_pro_list = []
all_cancer_cell_name_list = []
# key: cancer, value: (start_por, end_pot)
all_cancer_anno_start_end_pos_dict = {}
pos_count = -1
for cancer in cancer_cell_dict.keys():
    cell_name_list = cancer_cell_dict[cancer]
    start_pos = pos_count + 1
    for cell in cell_name_list:
        pos_count += 1
        all_cancer_cell_name_list.append(cell)
        self_pro = cell_selfValue_dict[cell]
        other_pro = 1 - self_pro
        self_pro_list.append(self_pro)
        other_pro_list.append(other_pro)
    end_pos = pos_count
    all_cancer_anno_start_end_pos_dict[cancer] = (start_pos, end_pos)

# Draw the stacked bar plot across all 18 cancers (Not used)
width = 0.7
fig, ax = plt.subplots(figsize = (50,7))
row_num = len(all_cancer_cell_name_list)
bottom = np.zeros(row_num)
p = ax.bar(all_cancer_cell_name_list, self_pro_list, width, label = "correctly predicted proportions", bottom = bottom, color = "red")
bottom += self_pro_list
p = ax.bar(all_cancer_cell_name_list, other_pro_list, width, label = "other proportions", bottom = bottom, color = "grey")
ax.legend(bbox_to_anchor=(1.01, 1.05))
plt.xticks(rotation = 90)
out_all_bar = in_path + "all_cancer_barplot.png"
print(ax.patches[-1].get_x())
plt.xlim(ax.patches[0].get_x()-0.3, ax.patches[-1].get_x()+1)
for cancer in all_cancer_anno_start_end_pos_dict.keys():
    trans = ax.get_xaxis_transform()
    (start_pos, end_pos) = all_cancer_anno_start_end_pos_dict[cancer]
    ax.annotate(cancer, xy=((start_pos+end_pos)/2, -0.45), xycoords=trans, ha="center", va="top")
    ax.plot([start_pos-0.3, end_pos+0.3],[-0.4,-0.4], color="k", transform=trans, clip_on=False)
plt.title("Cell line bulk RNA decomposition results across 18 cancers")
plt.ylabel("Proportions")
plt.tight_layout()
#plt.show()
plt.savefig(out_all_bar)
plt.close()

# output the correctly prediction proportions across 18 cancers
out_file = in_path + "all_cancer_self_proportions.txt"
fout = open(out_file, "w")
for cancer in cancer_cell_dict.keys():
    cell_name_list = cancer_cell_dict[cancer]
    for cell in cell_name_list:
        self_pro = cell_selfValue_dict[cell]
        fout.write("\t".join([cancer, cell, str(self_pro)]) + "\n")
fout.close()

# prepare data for boxplot
cancer_selfValue_dict = {}
for cancer in cancer_cell_dict.keys():
    cell_name_list = cancer_cell_dict[cancer]
    for cell in cell_name_list:
        self_pro = cell_selfValue_dict[cell]
        if cancer in cancer_selfValue_dict.keys():
            old_list = cancer_selfValue_dict[cancer]
            old_list.append(self_pro)
            cancer_selfValue_dict[cancer] = old_list
        else:
            cancer_selfValue_dict[cancer] = [self_pro]

### draw boxplot for the CCLE decomposed results of the 18 cancers
# data for the figure
data = []
cancer_order_list = []
for cancer in cancer_selfValue_dict.keys():
    cancer_order_list.append(cancer)
    one_cancer_data = cancer_selfValue_dict[cancer]
    data.append(one_cancer_data)

# get cancer list with descending order of mean proportions
cancer_data_dict_for_sort = {}
cancer_data_dict_for_sort["cancer_type"] = cancer_order_list
cancer_data_dict_for_sort["data_mean"] = list(map(np.mean, data))
df_for_sort = pd.DataFrame.from_dict(cancer_data_dict_for_sort)
df_for_sort.sort_values(by = ["data_mean"], ascending = False, inplace = True)
sorted_cancer_list = list(df_for_sort["cancer_type"])

# sorted data for the figure
sorted_data = []
for cancer in sorted_cancer_list:
    one_cancer_data = cancer_selfValue_dict[cancer]
    sorted_data.append(one_cancer_data)

xtick_pos_list = []
xtick_label_list = []
#for i in range(len(data)):
#    xtick_pos_list.append(i + 1)
#    cancer = cancer_order_list[i]
#    one_cancer_data = data[i]
#    data_num = len(one_cancer_data)
#    xtick_label_list.append(cancer + "(n = " + str(data_num) + ")")

for i in range(len(sorted_data)):
    xtick_pos_list.append(i + 1)
    cancer = sorted_cancer_list[i]
    one_cancer_data = sorted_data[i]
    data_num = len(one_cancer_data)
    xtick_label_list.append(cancer + "(n = " + str(data_num) + ")")
    
# draw the figure
all_fontsize = 14
in_title = "CCLE cell line decomposition results for 18 cancers"
out_file = in_path + "CCLE_cell_decomposed_proportions_for_18_cancers_boxplot.png"
plt.figure(figsize =(10, 8))
#plt.boxplot(data)
plt.boxplot(sorted_data)
plt.xticks(xtick_pos_list, xtick_label_list, fontsize = all_fontsize, rotation = 90)
plt.yticks(fontsize = all_fontsize)
plt.title(in_title, fontsize = all_fontsize)
plt.ylabel("Correctly predicted proportions", fontsize = all_fontsize)
plt.tight_layout()
plt.savefig(out_file)
plt.close()
