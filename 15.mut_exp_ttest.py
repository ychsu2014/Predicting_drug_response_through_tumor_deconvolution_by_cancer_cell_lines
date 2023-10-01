in_path = "/home/UTHSCSA_drug/data/mut_exp_data/"

in_cancer = "LungCancer"
in_RNA = in_path + in_cancer + "_RNA_overlapped_genes_geneSubset.h5ad"
#in_drug = in_path + in_cancer +"_addCancerType_addDrugName_None.txt"
in_idx = in_path + in_cancer + "_RNA_idx_sID.txt"
# afatinib: BRD-K66175015-001-09-0::2.5::HTS, alpelisib: BRD-K54997624-001-06-0::2.5::HTS
#in_drug_ID = "BRD-K66175015-001-09-0::2.5::HTS"


import anndata as ad
import pandas as pd
import numpy as np
import rpy2.robjects as robjects

def ttestByR(inData1, inData2):
	vector1 = robjects.FloatVector(inData1)
	vector2 = robjects.FloatVector(inData2)
	t_test = robjects.r["t.test"]
	result = t_test(vector1, vector2)
	p_value = float(np.array(result[2])[0])
	return(p_value)
