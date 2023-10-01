in_CCLE = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194.txt"
in_TCGA = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7781.txt"
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194_all_gene_num_by_cancer.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7781_all_gene_num_by_cancer.txt"


def calculateGeneNum(inFile, inGeneColnum, inCancerColnum, sep = "\t"):
    f = open(inFile)
    lines = f.readlines()
    lines = lines[1:]
    cancer_gene_dict = {}
    print(len(lines))
    count = 0
    for line in lines:
        count += 1
        if count % 10000 == 0:
            print(count)
        cols = line.strip("\n").split(sep)
        scRNA_cancer_type = cols[inCancerColnum]
        gene = cols[inGeneColnum]
        if scRNA_cancer_type in cancer_gene_dict.keys():
            gene_list = cancer_gene_dict[scRNA_cancer_type]
            gene_list.append(gene)
            gene_list = list(set(gene_list))
            cancer_gene_dict[scRNA_cancer_type] = gene_list
        else:
            cancer_gene_dict[scRNA_cancer_type] = [gene]
    return(cancer_gene_dict)

CCLE_cancer_gene_dict = calculateGeneNum(in_CCLE, 0, 32, ",")
TCGA_cancer_gene_dict = calculateGeneNum(in_TCGA, 6, 16)

fout1 = open(out_file1, "w")
fout2 = open(out_file2, "w")
for scRNA_cancer_type in CCLE_cancer_gene_dict.keys():
    print(scRNA_cancer_type)
    gene_list = CCLE_cancer_gene_dict[scRNA_cancer_type]
    gene_num = len(gene_list)
    fout1.write(scRNA_cancer_type + "\t" + str(gene_num) + "\n")

for scRNA_cancer_type in TCGA_cancer_gene_dict.keys():
    print(scRNA_cancer_type)
    gene_list = TCGA_cancer_gene_dict[scRNA_cancer_type]
    gene_num = len(gene_list)
    fout2.write(scRNA_cancer_type + "\t" + str(gene_num) + "\n")

fout1.close()
fout2.close()
