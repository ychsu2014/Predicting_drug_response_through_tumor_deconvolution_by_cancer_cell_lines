# script for volcano plot of all cancer types
in_file = "exp_lung19_all_cancers_adjusted_pvalues.txt"

out_file = "exp_all_cancer_volcano_plot.png"

in_gene_colnum = 1
in_drug_colnum = 3
in_X_colnum = 6
in_Y_colnum = 8
in_group_colnum = 0

import matplotlib
import matplotlib.pylab as plt
import random
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np

#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

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

# "Cancer_type\tgene\tprotein_hgvs\tdrug\tmut_num\twt_num\tmut_mean\twt_mean\tmut_wt_diff\tp_values\tadj_pvalues\n"
# "Cancer_type\tgene\tdrug\tmut_num\twt_num\tmut_mean\twt_mean\tmut_wt_diff\tp_value\tadj_p_value\n"
# inXColname: differences of the mean log fold change between mutant and wildtype
# inYColname: adjusted p-values
# inMidLabel: points at background
# inHighlightDict: keys: cancer type; values: labels for hightlight points
#def volcanoPlot(inDf, inXColname, inYColname, inGroupColname, outFile, inHighlightDict = cancer_type_dict):
def volcanoPlot(inXList, inYList, inGroupDict, outFile, inHighlightDict = cancer_type_dict, inColorDict = color_dict):
    plt.figure(figsize = (8, 80))  #8,80
    # s: marker size
    plt.scatter(x = inXList, y = inYList, s = 10)
    # generate a list of colors
    #groupNum = len(inHighlightDict.keys())
    ### get random colors from list of named colors
    #color_list = random.choices(list(mcolors.CSS4_COLORS.keys()), k = groupNum)
    ### other ways to generate random colors
    #color_list = ["#" + ''.join([random.choice('0123456789ABCDEF') for i in range(6)]) for j in range(groupNum)]
    # highlight different cancer types
    for i in range(len(inHighlightDict.keys())):
        cancer_type = list(inHighlightDict.keys())[i]
        #print(cancer_type)
        cancer_label = inHighlightDict[cancer_type]
        fig_color = inColorDict[cancer_type]
        #cancer_df = inDf[inDf[inGroupColname] == cancer_type]
        if cancer_label in inGroupDict.keys():
            print(cancer_label)
            print(fig_color)
            (cancerX, cancerY) = inGroupDict[cancer_label]
            plt.scatter(x = cancerX, y = cancerY, s = 10, label = cancer_label, color = fig_color)
    #texts=[]
    #for i,r in up.iterrows():
    #    texts.append(plt.text(x=r['diffmeanlogFC'],y=-np.log10(r['adj.P.Val']),s=i))
    # lw: line width
    #adjust_text(texts,arrowprops=dict(arrowstyle="->", color='black', lw=0.5))
    plt.xlabel("Differences of the log2 fold change between mutant and wildtype")
    plt.ylabel("-log(adjusted p-values)")
    #plt.axvline(-2,color="grey",linestyle="--")
    #plt.axvline(2,color="grey",linestyle="--")
    #plt.axhline(2,color="grey",linestyle="--")
    #plt.legend(loc = "best")
    plt.legend(loc = "lower center" ,bbox_to_anchor = (0.5, -1.26), ncol = 3)
    plt.tight_layout()
    plt.show()
    #plt.savefig(outFile)

# key: cancer_type, value: [x_list, y_list]
cancer_x_y_dict = {}
total_x = []
total_y = []
f = open(in_file)
#lines = f.readlines()
#lines = lines[1:]
count = 0
for line in f:
    count += 1
    if count == 1:
        continue
    cols = line.strip("\n").split("\t")
    gene = cols[in_gene_colnum]
    drug = cols[in_drug_colnum]
    exp_high_low_diff = float(cols[in_X_colnum])
    adj_p_value = -np.log10(float(cols[in_Y_colnum]))
    cancer_type = cols[in_group_colnum]
    if gene != "" and drug != "":
        total_x.append(exp_high_low_diff)
        total_y.append(adj_p_value)
        if cancer_type in cancer_x_y_dict.keys():
            (cancer_x_list, cancer_y_list) = cancer_x_y_dict[cancer_type]
            cancer_x_list.append(exp_high_low_diff)
            cancer_y_list.append(adj_p_value)
            cancer_x_y_dict[cancer_type] = (cancer_x_list, cancer_y_list)
        else:
            cancer_x_list = [exp_high_low_diff]
            cancer_y_list = [adj_p_value]
            cancer_x_y_dict[cancer_type] = (cancer_x_list, cancer_y_list)

        
volcanoPlot(total_x, total_y, cancer_x_y_dict, out_file)
