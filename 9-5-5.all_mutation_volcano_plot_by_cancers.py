# volcano plot by cancer types
in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/mut_grouped_drug_prediction_t_test_p_values_all_mutations_all_cancer_adjusted_pvalues.txt"

out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/all_mutations_lung19/volcano_plots_by_cancertypes/"

in_mut_colnum = 2
in_X_colnum = 8
in_Y_colnum = 10
in_group_colnum = 0

import matplotlib.pylab as plt
import numpy as np

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

color_dict = {}
#color_list = ["red", "orange", "limegreen", "turquoise", "gold", "blueviolet", "peru", "sienna", "olive", "magenta", "palevioletred", "tan", "lightgreen", "deeppink", "grey", "darkviolet", "darkgreen", "royalblue", "navy"]
color_dict["BileDuctCancer"] = "red"
color_dict["BladderCancer"] = "orange"
color_dict["BrainCancer"] = "limegreen"
color_dict["BreastCancer"] = "turquoise"
color_dict["ColonColorectalCancer"] = "gold"
color_dict["EndometrialUterineCancer"] = "blueviolet"
color_dict["EsophagealCancer"] = "peru"
color_dict["GastricCancer"] = "sienna"
color_dict["HeadandNeckCancer"] = "olive"
color_dict["KidneyCancer"] = "magenta"
color_dict["LiverCancer"] = "palevioletred"
color_dict["LungCancer"] = "tan"
color_dict["OvarianCancer"] = "lightgreen"
color_dict["PancreaticCancer"] = "deeppink"
color_dict["ProstateCancer"] = "grey"
color_dict["Sarcoma"] = "darkviolet"
color_dict["SkinCancer"] = "darkgreen"
color_dict["ThyroidCancer"] = "royalblue"

# inXColname: differences of the mean log fold change between mutant and wildtype
# inYColname: adjusted p-values
# inMidLabel: points at background
# inHighlightDict: keys: cancer type; values: labels for hightlight points
def volcanoPlot(inGroupDict, outPath, inHighlightDict = cancer_type_dict, inColorDict = color_dict):
    # highlight different cancer types
    for i in range(len(inHighlightDict.keys())):
        cancer_type = list(inHighlightDict.keys())[i]
        print(cancer_type)
        cancer_label = inHighlightDict[cancer_type]
        fig_color = inColorDict[cancer_type]
        if cancer_type in inGroupDict.keys():
            (cancerX, cancerY) = inGroupDict[cancer_type]
            plt.scatter(x = cancerX, y = cancerY, s = 10, label = cancer_label, color = fig_color)
        plt.xlabel("Differences of the log2 fold change between mutant and wildtype")
        plt.ylabel("-log(adjusted p-values)")
        plt.legend(loc = "lower center" ,bbox_to_anchor = (0.5, -1.26), ncol = 3)
        plt.tight_layout()
        plt.savefig(outPath + cancer_type.replace(" ", "_") + "_volcano_plot_removeBlank.png")
        plt.close()

# key: cancer_type, value: [x_list, y_list]
cancer_x_y_dict = {}
f = open(in_file)
#lines = f.readlines()
#lines = lines[1:]
count = 0
for line in f:
    count += 1
    if count == 1:
        continue
    cols = line.strip("\n").split("\t")
    mut = cols[in_mut_colnum]
    mut_wt_diff = float(cols[in_X_colnum])
    adj_p_value = -np.log10(float(cols[in_Y_colnum]))
    cancer_type = cols[in_group_colnum]
    if mut != "":
        if cancer_type in cancer_x_y_dict.keys():
            (cancer_x_list, cancer_y_list) = cancer_x_y_dict[cancer_type]
            cancer_x_list.append(mut_wt_diff)
            cancer_y_list.append(adj_p_value)
            cancer_x_y_dict[cancer_type] = (cancer_x_list, cancer_y_list)
        else:
            cancer_x_list = [mut_wt_diff]
            cancer_y_list = [adj_p_value]
            cancer_x_y_dict[cancer_type] = (cancer_x_list, cancer_y_list)

volcanoPlot(cancer_x_y_dict, out_path)
