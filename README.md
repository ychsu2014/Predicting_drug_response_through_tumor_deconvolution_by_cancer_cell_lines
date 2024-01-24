# [Predicting drug response through tumor deconvolution by cancer cell lines](https://github.com/ychsu2014/Predicting_drug_response_through_tumor_deconvolution_by_cancer_cell_lines)
## Summary
Large-scale cancer drug sensitivity data have become available for a collection of cancer cell lines, but only limited drug response data from patients are available. Bridging the gap of pharmacogenomics knowledge between in vitro and in vivo datasets remains challenging. In this study, we trained a deep-learning model, Scaden-CA, for deconvoluting tumor data into proportions of cancer type-specific cell lines. Then, a drug response prediction method was developed using the deconvoluted proportions and the drug sensitivity data from cell lines. The Scaden-CA model showed excellent performance in terms of concordance correlation coefficients (> 0.9 for model testing) and the correctly deconvoluted rate (> 70% across most cancers) for model validation using CCLE bulk RNA data. We applied the model to tumors in the TCGA dataset and examined associations between predicted cell viability and mutation status or gene expression levels for understanding underlying mechanisms of potential value for drug repurposing.

## Datasets used in the study
* Single cell RNA (scRNA) data of cell lines: http://singlecell.broadinstitute.org/single_cell/study/SCP542/pan-cancer-cell-line-heterogeneity
* TCGA data：https://xenabrowser.net/datapages/
* CCLE data: https://depmap.org/portal/ (2022Q2 version)
* PRISM drug screening dataset: https://depmap.org/repurposing/
* TCGA clinical data: downloaded by TCGAbiolinks (R/Bioconductor, https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)
* OncoKB database: https://www.oncokb.org/
* COSMIC Cancer Gene Census (CGC, v97) website: https://cancer.sanger.ac.uk/census

## Descriptions of the codes in this repository
### Dataset preprocessing
* 1-1 ~ 3-8-2: scRNA/TCGA/CCLE/PRISM datasets preprocessing
### Scaden-CA model training/testing/validation/application
* 4-1: Bulk RNA data simulation for Scaden-CA model training/testing
* 4-2-1_4-2-4: Scaden-CA model training/testing
* 4-3: draw figures to visualize the loss of the Scaden-CA model
* 4-4-1 ~ 4-4-5: Apply Scaden-CA model to CCLE RNA-Seq datasets
* 4-5-1 ~ 4-5-4: Apply Scaden-CA model to TCGA RNA-Seq datasets
### Scaden-CA model assessment
* 5-1 ~ 5-5: Estimate the optimal number of cell lines for Scaden-CA models
* 5-6: Heatmap for the sigle cell RNA data of the lung cancer cell lines
* 5-7-1 ~ 5-7-9: Simulation method comparison
### Combining oncoKB and COSMIC data in cell line selection
* 6-1 ~ 6-7: TCGA/CCLE mutation data preprocessing
* 7-1-1 ~ 7-1-13: Get overlapped mutatations between CCLE/TCGA mutation data and the records in oncoKB database
* 7-2-1 ~ 7-2-6: Create matrix to record the CGC cancer driver genes which CCLE/TCGA mutation occurred on
* 8-1 ~ 8-4: Reduce the number of cell lines for Scaden-CA models and assess the model performance
* 8-5 ~ 8-8-8: Apply the model trained with reduced cell line for lung cancer to TCGA/CCLE lung cancer data
### Drug response prediction for TCGA tumors
* 9-1-1 ~ 9-2-2: Drug reponse prediction for TCGA tumors by using the deconvolution results
* The "9-3-1_all_mutations" folder: Split TCGA tumors by mutation status and perform t-test to compare the drug responses between wild-type and mutants
* The "9-3-2_all_CGC_genes" folder: Split TCGA tumors by mutation status of CGC cancer driver genes and perform t-test to compare the drug responses between wild-type and mutants
* 9-4-1 ~ 9-4-4: p-value adjustment for the t-tests of Cancer-Gene-Drug combinations and visualization
* 9-5-1 ~ 9-5-6: p-value adjustment for the t-tests of Cancer-[gene, AA change]-Drug combinations and visualization
* 9-6-1 ~ 9-6-4: Split TCGA tumors by high/low expression of genes and perform t-tests to compare the drug responses between high/low groups
* 9-7: Use TCGA clinical data for the validation of drug response prediction



