# columns: sample, gene1_name, gene2_name, ...
# row: sample_name, 1 or 0, 1 or 0, ...
#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194.txt"
in_CGC = "/home/CBBI/hsuy1/projects/drugResponse/data/COSMIC/Census_allMon_Jan_30_21_08_13_2023.tsv"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered_194_CGC_matrix.txt"

in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194_filteredLungCancer19.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/CCLE_all_mutations_filtered194_filteredLungCancer19_CGC_matrix.txt"

f = open(in_file)
f_gene = open(in_CGC)
fout = open(out_file, "w")

lines = f_gene.readlines()
lines = lines[1:]
CGC_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    gene = cols[0]
    CGC_list.append(gene)

lines = f.readlines()
header = lines[0]
lines = lines[1:]
sID_mutGene_dict = {}
for line in lines:
    cols = line.strip("\n").split(",")
    sID = cols[15]
    gene = cols[0]
    if sID in sID_mutGene_dict.keys():
        old_list = sID_mutGene_dict[sID]
        old_list.append(gene)
        old_list = list(set(old_list))
        sID_mutGene_dict[sID] = old_list
    else:
        sID_mutGene_dict[sID] = [gene]

# write header
fout.write("cell_line")
for gene in CGC_list:
    fout.write("\t" + gene)
fout.write("\n")

for sID in sID_mutGene_dict.keys():
    gene_list = sID_mutGene_dict[sID]
    fout.write(sID)
    for gene in CGC_list:
        if gene in gene_list:
            fout.write("\t1")
        else:
            fout.write("\t0")
    fout.write("\n")

f.close()
f_gene.close()
fout.close()
