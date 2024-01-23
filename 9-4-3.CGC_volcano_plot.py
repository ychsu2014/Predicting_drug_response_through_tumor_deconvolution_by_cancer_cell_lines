# dataframe with:
# (1) differences of predicted log2 fold change between mutant and wildtype
# (2) p-values
# (3) adjusted p-values 
# (4) cancer types
# (5) mutation
# (6) drug

#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/mut_grouped_drug_prediction_t_test_p_values_all_in_oncoKB_v2_adjusted.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/cancer_mutation_drug_plots/mut_grouped_drug_prediction_t_test_p_values_all_in_oncoKB_v2_adjusted.png"
#in_X_colname = "mut_wt_diff"
#in_Y_colname = "adj_pvalues"
#in_group_colname = "Cancer_type"

in_file = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/drug_score/exponential_log/mut_grouped_drug_prediction_t_test_p_values_CGC_all_cancer_adjusted_pvalues.txt"
#"/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/mut_grouped_drug_prediction_t_test_p_values_CGC_all_cancer_adjusted_pvalues.txt"
out_file = "test.png"
"/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/drug_score/exponential_log/mut_grouped_drug_prediction_t_test_p_values_CGC_all_cancer_adjusted_pvalues_v2.png"
in_X_colname = "mut_wt_diff"
in_Y_colname = "adj_p_value"
in_group_colname = "Cancer_type"

import matplotlib.pylab as plt
import random
import matplotlib.colors as mcolors
import pandas as pd
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

# "Cancer_type\tgene\tprotein_hgvs\tdrug\tmut_num\twt_num\tmut_mean\twt_mean\tmut_wt_diff\tp_values\tadj_pvalues\n"
# "Cancer_type\tgene\tdrug\tmut_num\twt_num\tmut_mean\twt_mean\tmut_wt_diff\tp_value\tadj_p_value\n"
# inXColname: differences of the mean log fold change between mutant and wildtype
# inYColname: adjusted p-values
# inMidLabel: points at background
# inHighlightDict: keys: cancer type; values: labels for hightlight points
def volcanoPlot(inDf, inXColname, inYColname, inGroupColname, outFile, inHighlightDict = cancer_type_dict):
    plt.figure(figsize = (8, 80))
    # s: marker size
    plt.scatter(x = inDf[inXColname], y = inDf[inYColname].apply(lambda x:-np.log10(x)), s = 10)
    # generate a list of colors
    #groupNum = len(inHighlightDict.keys())
    color_list = ["red", "orange", "limegreen", "turquoise", "gold", "blueviolet", "peru", "sienna", "olive", "magenta", "palevioletred", "tan", "lightgreen", "deeppink", "grey", "darkviolet", "darkgreen", "royalblue", "navy"]
    ### get random colors from list of named colors
    #color_list = random.choices(list(mcolors.CSS4_COLORS.keys()), k = groupNum)
    ### other ways to generate random colors
    #color_list = ["#" + ''.join([random.choice('0123456789ABCDEF') for i in range(6)]) for j in range(groupNum)]
    # highlight different cancer types
    for i in range(len(inHighlightDict.keys())):
        cancer_type = list(inHighlightDict.keys())[i]
        print(cancer_type)
        cancer_label = inHighlightDict[cancer_type]
        fig_color = color_list[i]
        cancer_df = inDf[inDf[inGroupColname] == cancer_type]
        plt.scatter(x = cancer_df[inXColname], y = cancer_df[inYColname].apply(lambda x:-np.log10(x)), s = 10, label = cancer_label, color = fig_color)
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
    plt.savefig(outFile)

# for testing
#test = {}
#test["Cancer_type"] = ["PancreaticCancer", "Sarcoma", "SkinCancer", "SkinCancer", "SkinCancer"]
#test["adj_pvalues"] = [0.001, 0.02, 0.00001, 0.0000006, 0.9]
#test["fold_change"] = [-0.2, 0.3, -0.1, 0.07, 0.8]
#test_df = pd.DataFrame.from_dict(test)
#volcanoPlot(test_df, "fold_change", "adj_pvalues", "Cancer_type", "test.png")

df = pd.read_csv(in_file, sep = "\t")
volcanoPlot(df, in_X_colname, in_Y_colname, in_group_colname, out_file)
