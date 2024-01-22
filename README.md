# [Predicting drug response through tumor deconvolution by cancer cell lines](https://github.com/ychsu2014/Predicting_drug_response_through_tumor_deconvolution_by_cancer_cell_lines)
## Summary
Large-scale cancer drug sensitivity data have become available for a collection of cancer cell lines, but only limited drug response data from patients are available. Bridging the gap of pharmacogenomics knowledge between in vitro and in vivo datasets remains challenging. In this study, we trained a deep-learning model, Scaden-CA, for deconvoluting tumor data into proportions of cancer type-specific cell lines. Then, a drug response prediction method was developed using the deconvoluted proportions and the drug sensitivity data from cell lines. The Scaden-CA model showed excellent performance in terms of concordance correlation coefficients (> 0.9 for model testing) and the correctly deconvoluted rate (> 70% across most cancers) for model validation using CCLE bulk RNA data. We applied the model to tumors in the TCGA dataset and examined associations between predicted cell viability and mutation status or gene expression levels for understanding underlying mechanisms of potential value for drug repurposing.

## Datasets used in the study
* Single cell RNA (scRNA) data of cell lines: http://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity \
* TCGA data：https://xenabrowser.net/datapages/ \
* CCLE data: https://depmap.org/portal/ (2022Q2 version) \
* PRISM drug screening dataset: https://depmap.org/repurposing/

## Descriptions of the codes in this repository
* Dataset preprocessing
  * 1-1 ~ 3-8-1: scRNA/TCGA/CCLE/PRISM datasets preprocessing\
* Scaden-CA model training/testing/validation/application
  * 4-1: Bulk RNA data simulation for Scaden-CA model training/testing\
  * 4-2-1_4-2-4: Scaden-CA model training/testing\
  * 4-3: draw model loss figure\
  * 4-4-1 ~ 4-7-2: Apply Scaden-CA model to CCLE and TCGA RNA-Seq datasets\
* Incoporating TCGA/CCLE mutation data for cell line reduction
  * 5-1 ~ 5-7: TCGA/CCLE mutation data preprocessing\
  * 6-1 ~ 6-7: Get overlapped mutatations between CCLE/TCGA mutation data and the records in oncoKB database\
  * 7-1 ~ 7-6: Create matrix to record the CGC cancer driver genes which CCLE/TCGA mutation occurred on\
  * 8-1 ~ 8-4: Reduce the number of cell lines for Scaden-CA models and assess the model performance
  * 8-5 ~ 8-8-8: Apply the model trained with reduced cell line for lung cancer to TCGA/CCLE lung cancer data




