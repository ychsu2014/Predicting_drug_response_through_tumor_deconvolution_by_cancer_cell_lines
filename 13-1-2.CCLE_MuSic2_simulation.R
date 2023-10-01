#argvs <- commandArgs(trailingOnly = TRUE)

# ex. Breast_Cancer
#cancer_type <- argvs[2]

rm(list = ls())

cancer_type <- "Thyroid_Cancer"
#cancer_type <- "Skin_Cancer"
#cancer_type <- "Sarcoma" # std n-1 can not run
#cancer_type <- "Prostate_Cancer" # std n-1 can not run
#cancer_type <- "Pancreatic_Cancer"
#cancer_type <- "Ovarian_Cancer"
#cancer_type <- "Lung_Cancer"
#cancer_type <- "Liver_Cancer" # std n-1 can not run
#cancer_type <- "Kidney_Cancer"
#cancer_type <- "Head_and_Neck_Cancer"
#cancer_type <- "Gastric_Cancer" # std n-1 can not run
#cancer_type <- "Esophageal_Cancer"
#cancer_type <- "Endometrial_Uterine_Cancer"
#cancer_type <- "Colon_Colorectal_Cancer"
#cancer_type <- "Brain_Cancer"
#cancer_type <- "Bladder_Cancer"
#cancer_type <- "Bile_Duct_Cancer" # std n-1 can not run
#cancer_type <- "Breast_Cancer"

cancer_type2 <- gsub("_", "", cancer_type)

in_data <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/keep_only_TCGA_genes_lung_19/", cancer_type, "_norm_counts_all.txt", sep = "")

in_phe <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/scRNA/preprocessing_onlyNormalized/data_for_MuSiC/", cancer_type, "_cellLine_cancertypes_sID.txt", sep = "")

in_prop <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/", cancer_type2, "_RNA_ovelapped_genes_geneSubset_prediction.txt", sep = "")

in_cellLine_idx <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/idx_files/", cancer_type2, "_RNA_idx_sID.txt", sep = "")

out_file1 <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/", cancer_type, "_sampled_proportion.txt", sep = "")

out_file2 <- paste("/Users/ychsu/Library/CloudStorage/GoogleDrive-ychsu20130517@gmail.com/Other computers/My Laptop/TIGP_bioinformatics/research/UTHSCSA/drugResponse/data/ICIBM_2023/CCLE/RNA/by_cancer/h5ad/pseudo_bulk_RNA/MuSiC2_simulation/", cancer_type, "_pseudo_bulk.txt", sep = "")

exp <- as.matrix(read.table(in_data, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE))
exp <- t(exp)
phe <- read.table(in_phe, row.names = 1, header = TRUE, sep = "\t")
prop <- read.table(in_prop, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
prop_df <- data.frame(prop, check.names = FALSE)
prop_df <- cbind(cellID = rownames(prop_df), prop_df)
cell_line_idx <- read.table(in_cellLine_idx, header = FALSE, sep = "\t")
colnames(cell_line_idx) <- c("cellID", "cell_name", "cell_name2", "cell_name3", "depmap_ID")
prop_cell_df <- merge(prop_df, cell_line_idx, by = "cellID")

# use 500 cells for one bulk RNA sample
cell_num <- 500
metadata_name <- paste(cancer_type, "_SCE", sep = "")

library("SingleCellExperiment")
library("MCMCpack")
library("stats")
library("dampack")
library("GMCM")
#library("rBeta2009")
#library("dplyr")
#require(MASS)

gene_metadata <- data.frame(
  geneName = rownames(exp)
)
cell_metadata <- data.frame(
  sampleID = phe["sampleID"],
  cellLine = phe["cellLine"],
  cancerType = phe["cancertype"]
)
exp_SCE <- SingleCellExperiment(
  assay = list(counts = exp),
  colData = cell_metadata,
  rowData = gene_metadata,
  metadata = list(name = metadata_name)
)

### get the subject level variation
#music_Sigma = function(x, non.zero, markers, clusters, samples, select.ct){
x <- exp_SCE
non.zero <- FALSE
markers <- rownames(exp_SCE@assays@data$counts)
clusters <- "cellLine"
# sampleID is the ID of each cell lines
samples <- "sampleID"
select.ct <- unique(exp_SCE@colData$cellLine)

if(!is.null(select.ct)){
  x = x[, x@colData[, clusters] %in% select.ct]
}
if(non.zero){  ## eliminate non expressed genes
  x <- x[rowSums(counts(x))>0, ]
}
  
clusters <- as.character(colData(x)[, clusters])
samples <- as.character(colData(x)[, samples])

Sigma <- sapply(unique(clusters), function(ct){
  apply(sapply(unique(samples), function(sid){
    y = counts(x)[,clusters %in% ct & samples %in% sid]
    if(is.null(dim(y))){
      temp_values1 <- y/sum(y)
      temp_values1 <- replace(temp_values1, is.na(temp_values1), 0)
      return(temp_values1)
    }else{
      temp_values2 <- rowSums(y)/sum(y)
      temp_values2 <- replace(temp_values2, is.na(temp_values2), 0)
      return(temp_values2)
    }
  }), 1, var, na.rm = TRUE)
})
  
if(!is.null(select.ct)){
  m.ct = match(select.ct, colnames(Sigma))
  Sigma = Sigma[, m.ct]
}
  
if (!is.null(markers)){
  ids <- intersect(unlist(markers), rownames(x))
  m.ids = match(ids, rownames(x))
  Sigma <- Sigma[m.ids, ]
}
#  return(Sigma = Sigma)
#}
#sub.var matrix, 2 by G, scale of log-normal distribution for each subjects at genes. Part of them are zeros.
sub.var <- t(Sigma)

### create mu: matrix, 2 by G, the mean expression of genes in each cell types.
scRNA_count <- data.frame(exp_SCE@assays@data$counts, check.names = FALSE)
scRNA_count <- t(scRNA_count)
scRNA_count2 <- cbind(cellID = rownames(scRNA_count), scRNA_count)

cellID <- rownames(exp_SCE@colData)
sampleID <- exp_SCE@colData$sampleID
cellLine <- exp_SCE@colData$cellLine
cancertype <- exp_SCE@colData$cancertype
phenotype_data <- data.frame(cellID, sampleID, cellLine, cancertype)
mergedData <- merge(phenotype_data, scRNA_count2, by = "cellID")

cell_line_list = unique(cellLine)
count_num <- 0
for (one_cell_line in cell_line_list){
  count_num <- count_num + 1
  temp_data <- mergedData[mergedData["cellLine"] == one_cell_line,]
  temp_data <- temp_data[,5:length(colnames(temp_data))]
  temp_data <- data.frame(sapply(temp_data, function(x) as.numeric(as.character(x))))
  temp_data <- colMeans(temp_data)
  #####
  if (one_cell_line == "42MGBA"){
    # for brain cancer
    one_cell_line = "X42MGBA"
  }else if (one_cell_line == "8305C"){
    # for thyroid cancer
    one_cell_line = "X8305C"
  }else if (one_cell_line == "2313287"){
    # for gastric cancer
    one_cell_line = "X2313287"
  }
  #####
  if (count_num == 1){
    # R variable name can not start with numbers?
    eval(parse(text = paste("mu <- data.frame(", one_cell_line, " = temp_data, check.names = FALSE)")))
  }else{
    eval(parse(text = paste("mu <- data.frame(mu, ", one_cell_line," = temp_data, check.names = FALSE)")))  
  }
}
mu <- t(mu)
#####
mu_rownames <- rownames(mu)
new_mu_rownames <- c()
for (i in rownames(mu)){
  if (i == "X42MGBA"){
    # for brain cancer
    i = "42MGBA"
  }else if (i == "X8305C"){
    # for thyroid cancer
    i = "8305C"
  }else if (i == "X2313287"){
    # for gastric cancer
    i = "2313287"
  }
  new_mu_rownames <- c(new_mu_rownames, i)
}
row.names(mu) <- new_mu_rownames
#####
### estimate the parameter for the dirichlet distribution of the cell type proportion
# pro_mena and pro_std are from the deconvoluted cell line proportions of breast cancers
#pro_celltypes <- c("MDAMB436", "T47D", "BT474", "BT549", "HMC18", "EFM192A", "HCC1428", "HCC1419", "ZR751", "HDQP1", "CAMA1", "MCF7", "MDAMB361", "HCC38")
end_num <- 1+length(cell_line_list)
pro_mean <- colMeans(prop_cell_df[,2:end_num])
# in R, the std is calculated with denominator with n-1 (sample), but in estimating the dirichlet parameters, use denomimator with n (population)
pro_std <- colSds(data.matrix(prop_cell_df[, 2:end_num])) *sqrt((length(cell_line_list)-1)/length(cell_line_list))
print(pro_std)
#pro_mean <- c(0.07515139162857144, 0.09224497230714286, 0.05651059504285715, 0.07793501062142857, 0.07458123421428571, 0.06092830560714285, 0.051675554985714287, 0.06609613957857143, 0.03169527396428572, 0.06706671692142857, 0.1334181267142857, 0.06833162639285714, 0.07241234108571429, 0.07195269159285714)
#pro_std <- c(0.2376770283425068, 0.2310184837037112, 0.1525124176875083, 0.24078526804595246, 0.22612998901380216, 0.1675004109160061, 0.1499176445289306, 0.17589121847835462, 0.08443575323306221, 0.21326577350573486, 0.23750004924830356, 0.18180897346435634, 0.1715405780582676, 0.22737824232534518)
diri_para <- dirichlet_params(pro_mean, pro_std)

### simulation
n.sc <- length(cellID)
n.bulk <- length(cell_line_list)
G <- length(markers)
sigma <- NULL
sort.bulk <- FALSE

# generate cell number and cell type proportion
K <- length(cell_line_list)
N <- n.sc + n.bulk
Ni = rep(500, N)
set.seed(100)
p = rdirichlet(N, diri_para)
Nk = round(p[,1] * Ni)
for (i in 2:n.bulk){
  Nk <- cbind(Nk, round(p[,i] * Ni))
}
p <- Nk/Ni
# sort p and Nk
colnames(p) = colnames(prop_cell_df[,2:end_num])
colnames(Nk) = colnames(prop_cell_df[,2:end_num])
p <- p[,rownames(mu)]
Nk <- Nk[,rownames(mu)]

#p_df <- data.frame(p)
#colnames(p_df) <- colnames(prop_cell_df[,2:15])
  
# Generate subject level variation
count_num <- 0
for (j in 1:length(rownames(sub.var))){
  count_num <- count_num + 1
  print(count_num)
  if (count_num == 1){
    sub.Fac <- list(sapply(1:G, function(i){rlnorm(N, 0, sub.var[j, i])}))
  }else{
    temp_data <- list(sapply(1:G, function(i){rlnorm(N, 0, sub.var[j, i])}))
    sub.Fac <- c(sub.Fac, temp_data)
  }
}

Bulk.null = FALSE; SC.null = FALSE;
if(n.bulk == 0) Bulk.null = TRUE;
if(n.sc == 0) SC.null = TRUE;

if(is.null(sigma)){
  if(!Bulk.null){
    Xjg = NULL    # generate bulk tissue data
    # i: number of bulk samples
    for(i in 1:n.bulk){
      RNA_data_temp <- rep(0, length(markers))
      # j: number of cell lines
      for(j in 1:n.bulk){
        # there are zeros in Nk, skip zeros for bulk tissue data generation
        if (Nk[i,j] > 1){
          samp_num_temp <- Nk[i,j]
          #print(samp_num_temp)
          RNA_data_temp <- RNA_data_temp + colSums(sapply(mu[j,] * sub.Fac[[j]][i, ], rpois, n = samp_num_temp))
        }
      }
      #print(RNA_data_temp)
      Xjg = cbind(Xjg, RNA_data_temp)
    }
  }
}

if(!Bulk.null){
  p.bulk = p[1:n.bulk, ];
  #if(sort.bulk){
  #  op.bulk = order(p.bulk[,1])
  #  p.bulk = p.bulk[op.bulk, ]; rownames(p.bulk) = paste0('Bulk', 1:n.bulk)
  #  Xjg = Xjg[, op.bulk]; colnames(Xjg) <- paste0('Bulk', 1:n.bulk)
  #}else{
  rownames(p.bulk) <- rownames(mu)
  colnames(Xjg) <- rownames(mu)
  #}
  #rownames(Xjg) <- GN;
  bulk.mtx = data.matrix(Xjg)
}

write.table(p.bulk, file = out_file1, sep = "\t")
write.table(bulk.mtx, file = out_file2, sep = "\t")
