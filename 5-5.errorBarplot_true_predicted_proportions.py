# output mean error barplots with standard error as error bar (all cancer merged and by cancer types)
in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/"
lung_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/lung_cancer_19/"
#####
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/figures_error_barplots/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/figures_error_barplots_interval0.03/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/figure_error_barplots_stderror/interval_0.1_0_1/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/error_proportion_estimation/simulation_data/figure_error_barplots_stderror/interval_0.01_0_0.1/"
out_file = out_path + "error_mean_stderror_all_cancers.txt"
out_file2 = out_path + "error_mean_stderror_all_cancers_combined.txt"

# tried different conditions: 0-1, 0-0.1
fig_lower_limit = 0.1
fig_upper_limit = 1.05
fig_step = 0.1
round_digit_num = 1
sample_num = 40
#fig_lower_limit = 0.01
#fig_upper_limit = 0.105
#fig_step = 0.01
#round_digit_num = 2
#sample_num = 30

import os
import pandas as pd
import anndata as ad
import math
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import stats

all_files = os.listdir(in_path)
pred_file_list = []
for afile in all_files:
    if "prediction" in afile:
        pred_file_list.append(afile)

upper_interval_dict = {}
##### try different intervals
#for i in np.arange(0.1, 1.05, 0.1).round(1):
#    upper_interval_dict[i] = str(round(i-0.1,1)) + "-" + str(round(i,1))
#for i in np.arange(0.01, 0.205, 0.01).round(2):
#    upper_interval_dict[i] = str(round(i-0.01, 2)) + "-" + str(round(i,2))
for i in np.arange(fig_lower_limit, fig_upper_limit, fig_step).round(round_digit_num):
    upper_interval_dict[i] = str(round(i-fig_step, round_digit_num)) + "-" + str(round(i, round_digit_num))

# inUpper: the upper limit of the proportion interval (0.1 -> 0-0.1, 0.2 -> 0.1-0.2, 0.3 -> 0.2-0.3, ...)
# inGroup: "T" -> True proportion, "P" -> predicted proportion
# inValue: the proportion that needs to be added to the dict
def addToDict(inDict, inUpper, inGroup, inValue):
    inKey = (inUpper, inGroup)
    if inKey in inDict.keys():
        oldList = inDict[inKey]
        oldList.append(inValue)
        inDict[inKey] = oldList
    else:
        inDict[inKey] = [inValue]
    return(inDict)

def addBothValueToDict(inDict, inTrueValue, inPredValue):
    # determine the interval and the upper limit
    #####
    inUpper = math.ceil(inTrueValue * int(1/fig_step)) / (1/fig_step)
    #inUpper = math.ceil(inTrueValue * 100) / 100.0
    # add true proportion to the dict
    inDict = addToDict(inDict, inUpper, "T", inTrueValue)
    #  add predicted proportion to the dict
    inDict = addToDict(inDict, inUpper, "P", inPredValue)
    return(inDict)

# key: uppper limit, value: list of errors
all_cancer_error_dict = {}
fout = open(out_file, "w")
fout.write("cancer_type\tinterval_upper\terror_mean\terror_se\terror_values\n")
for pred_file in pred_file_list:
    cancer_type = pred_file[:-15]
    #simu_file = "testing_" + cancer_type + "_processed.h5ad"
    simu_file = cancer_type + "_processed.h5ad"
    if cancer_type == "Lung_Cancer":
        pred_df = pd.read_csv(lung_path + pred_file, sep = "\t", index_col = 0)
        true_data = ad.read_h5ad(lung_path + simu_file)
        true_df = true_data.obs
    else:
        pred_df = pd.read_csv(in_path + pred_file, sep = "\t", index_col = 0)
        true_data = ad.read_h5ad(in_path + simu_file)
        true_df = true_data.obs
    # key: (range upper limit of the true proportion, "T"/"P"), value: [proportion]
    # T: true, P: predicted
    cancer_temp_dict = {}
    print(cancer_type)
    for cell_line in pred_df.columns:
        for s_idx in pred_df.index:
            pred_value = pred_df[cell_line].loc[s_idx]
            true_value = true_df[cell_line].loc[str(s_idx)]
            #print(cell_line)
            #print(s_idx)
            #print(pred_value)
            #print(true_value)
            cancer_temp_dict = addBothValueToDict(cancer_temp_dict, true_value, pred_value)
    interval_list = []
    interval_mean_list = []
    interval_se_list = []
    #####
    #for i in np.arange(0.1, 1.05, 0.1).round(1):
    #for i in np.arange(0.01, 0.205, 0.01).round(2):
    #for i in np.arange(0.1, 1.005, 0.01).round(2):
    for i in np.arange(fig_lower_limit, fig_upper_limit, fig_step).round(round_digit_num):
        true_value_list = cancer_temp_dict[(i, "T")]
        pred_value_list = cancer_temp_dict[(i, "P")]
        ##### select only 40 values for true/prediction group
        value_idx = random.sample(range(len(true_value_list)), k = sample_num)
        true_value_list = list(np.array(true_value_list)[value_idx])
        pred_value_list = list(np.array(pred_value_list)[value_idx])
        #####
        if len(true_value_list) == len(pred_value_list):
            pair_num = len(true_value_list)
        else:
            print("The number of true/predicted values is not consistent.")
        ### data for all cancers
        if i in all_cancer_error_dict.keys():
            all_cancer_error_list_temp  = all_cancer_error_dict[i]
        else:
            all_cancer_error_list_temp = []
        ###
        error_list = []
        for j in range(len(true_value_list)):
            true_temp = true_value_list[j]
            pred_temp = pred_value_list[j]
            error = abs((pred_temp-true_temp)/true_temp)
            #error = abs((pred_temp-true_temp))
            error_list.append(error)
            all_cancer_error_list_temp.append(error)
        ### data for all cancers
        all_cancer_error_dict[i] = all_cancer_error_list_temp
        ###
        # data for individual cancer types
        error_mean = np.mean(error_list)
        error_se = stats.sem(error_list)
        interval_temp = upper_interval_dict[i]
        interval_list.append(interval_temp)
        interval_mean_list.append(error_mean)
        interval_se_list.append(error_se)
        fout.write(cancer_type + "\t" + str(i) + "\t" + str(error_mean) + "\t" + str(error_se) + "\t" + str(error_list) + "\n")
    x_pos = np.arange(len(interval_list))
    plt.bar(x_pos, interval_mean_list, yerr = interval_se_list, capsize = 10)
    plt.ylabel("Error", fontsize = 14)
    plt.xlabel("Intervals of true proportions", fontsize = 14)
    plt.title(cancer_type, fontsize = 14)
    plt.xticks(x_pos, interval_list, rotation = 90)
    plt.tight_layout()
    plt.savefig(out_path + cancer_type + "_error_barplot.pdf", format = "pdf")
    plt.close()
fout.close()

### figure for all cancers
all_cancer_interval_list = []
all_cancer_interval_mean_list = []
all_cancer_interval_se_list = []
fout2 = open(out_file2, "w")
fout2.write("interval_upper\terror_mean\terror_se\terror_values\n")
for i in np.arange(fig_lower_limit, fig_upper_limit, fig_step).round(round_digit_num):
    interval_temp = upper_interval_dict[i]
    all_cancer_error_list = all_cancer_error_dict[i]
    all_cancer_error_mean = np.mean(all_cancer_error_list)
    all_cancer_error_se = stats.sem(all_cancer_error_list)
    all_cancer_interval_list.append(interval_temp)
    all_cancer_interval_mean_list.append(all_cancer_error_mean)
    all_cancer_interval_se_list.append(all_cancer_error_se)
    fout2.write("\t".join([str(i), str(all_cancer_error_mean), str(all_cancer_error_se), str(all_cancer_error_list)]) + "\n")
fout2.close()
x_pos = np.arange(len(all_cancer_interval_list))
plt.bar(x_pos, all_cancer_interval_mean_list, yerr = all_cancer_interval_se_list, capsize = 10)
plt.ylabel("Error", fontsize = 14)
plt.xlabel("Intervals of true proportions", fontsize = 14)
plt.title("All 18 cancer types", fontsize = 14)
plt.xticks(x_pos, all_cancer_interval_list, rotation = 90)
plt.tight_layout()
plt.savefig(out_path + "all_cancers_error_barplot.pdf", format = "pdf")
plt.close()
