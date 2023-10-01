in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/TCGA_CCLE_in_oncoKB_uni_var.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/TCGA_CCLE_in_oncoKB_uni_var_geneCount.txt"

file_colnum = 0
gene_colnum = 6

f = open(in_file)
fout = open(out_file, "w")

lines = f.readlines()
# key: filename, value: list of genes
file_gene_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    filename = cols[file_colnum]
    gene = cols[gene_colnum]
    if filename in file_gene_dict.keys():
        old_list = file_gene_dict[filename]
        old_list.append(gene)
        old_list = list(set(old_list))
        file_gene_dict[filename] = old_list
    else:
        file_gene_dict[filename] = [gene]

f.close()

for filename in file_gene_dict.keys():
    gene_list = file_gene_dict[filename]
    fout.write(filename + "\t" + str(len(gene_list)) + "\t" + "\t".join(gene_list) + "\n")

fout.close()
