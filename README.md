# [Predicting drug response through tumor deconvolution by cancer cell lines](https://github.com/ychsu2014/Predicting_drug_response_through_tumor_deconvolution_by_cancer_cell_lines)
## Summary
Large-scale cancer drug sensitivity data have become available for a collection of cancer cell lines, but only limited drug response data from patients are available. Bridging the gap of pharmacogenomics knowledge between in vitro and in vivo datasets remains challenging. In this study, we trained a deep-learning model, Scaden-CA, for deconvoluting tumor data into proportions of cancer type-specific cell lines. Then, a drug response prediction method was developed using the deconvoluted proportions and the drug sensitivity data from cell lines. The Scaden-CA model showed excellent performance in terms of concordance correlation coefficients (> 0.9 for model testing) and the correctly deconvoluted rate (> 70% across most cancers) for model validation using CCLE bulk RNA data. We applied the model to tumors in the TCGA dataset and examined associations between predicted cell viability and mutation status or gene expression levels for understanding underlying mechanisms of potential value for drug repurposing.

## Datasets used in the study
Single cell RNA (scRNA) data of cell lines: http://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity\
The Cancer Genome Atlas (TCGA) RNA, somatic mutation, and phenotypic data：https://xenabrowser.net/datapages/\
Cancer Cell Line Encyclopedia (CCLE) RNA data: https://depmap.org/portal/ (2022Q2 version)\
Profiling Relative Inhibition Simultaneously in Mixtures (PRISM) drug screening dataset: https://depmap.org/repurposing/\

## Codes description
1-1 ~ 2-5: scRNA data preprocessing\
3-1 ~ 3-4, 3-8-2: TCGA data preprocessing\
3-5, 3-6-2: CCLE data preprocessing\
3-6-1: scRNA/CCLE/drug datasets preprocessing\
3-7-1, 3-7-2: drug dataset preprocessing\
3-8-1: CCLE/TCGA data preprocessing\
4-1: Bulk RNA data simulation for Scaden-CA model training/testing\
4-2-1_4-2-4: Scaden-CA model training/testing\
4-3: draw model loss figure\
4-4-1, 4-4-2: CCLE/TCGA data preprocessing\
4-5-1, 4-5-2: Apply the Scaden-CA model to TCGA tumors\



