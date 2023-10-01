# this script is to reformat the data to the format that can be loaded using the "sc.read_10x_mtx" command
in_file = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_data.txt"
# gene list from the PBMC3k dataset used in the scanpy tutorial
in_gene = "/home/yuching/projects/drugResponse/data/scRNA/PBMC3K_geneID/features.tsv"

out_file1 = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_files/matrix.mtx"
out_file2 = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_files/barcodes.tsv"
out_file3 = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_files/features.tsv"

#from pyensembl import EnsemblRelease

fout_mat = open(out_file1, "w")
fout_barc = open(out_file2, "w")
fout_feat = open(out_file3, "w")
f_gene = open(in_gene)
f_data = open(in_file)

# geneName_ID_dict: key: gene_name, value: gene_id
lines = f_gene.readlines()
geneName_ID_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    gene_id = cols[0]
    gene_name = cols[1]
    geneName_ID_dict[gene_name] = gene_id
f_gene.close()
print("gene parsed.")

# first line: "GENE\tID1\tID2\t...\n"
lines = f_data.readlines()
print("data loaded.")
header = lines[0].strip("\n").split("\t")
sID_list = header[1:]
lines = lines[1:]
gene_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    gene_list.append(cols[0])

for sID in sID_list:
    fout_barc.write(sID + "\n")
fout_barc.close()
print("barcode done.")
# import Ensembl human release 107
#Ens107 = EnsemblRelease(107)
# write to "features.tsv"
for gene in gene_list:
    # get gene id by gene names
    #gene_id = Ens107.genes_by_name(gene)[0].gene_id
    gene_id = geneName_ID_dict[gene]
    fout_feat.write(gene_id + "\t" + gene + "\tGene Expression\n")
fout_feat.close()
print("feature done.")
# write to "matrix.mtx"
# header
fout_mat.write("%%MatrixMarket matrix coordinate real general\n")
fout_mat.write("%\n")
#nonzero_entry_num = len(np.where(data != 0)[0])
#fout_mat.write(str(len(gene_list)) + " " + str(len(sID_list)) + " " + str(nonzero_entry_num) + "\n")
row_count = 0
non_empty_count = 0
# (gene_idx, sID_idx)   exp_value
entry_dict = {}
for line in lines:
    row_count += 1
    cols = line.strip("\n").split("\t")
    gene_name = cols[0]
    exp_values = cols[1:]
    col_count = 0
    for exp_value in exp_values:
        col_count += 1
        if exp_value != "0":
            non_empty_count += 1
            entry_dict[(row_count, col_count)] = exp_value
print("parse to dict done.")

fout_mat.write(str(len(gene_list)) + " " + str(len(sID_list)) + " " + str(non_empty_count) + "\n")

p_count = 0
for key in entry_dict.keys():
    p_count += 1
    if p_count % 10000 ==0:
        print(p_count)
    row_c = key[0]
    col_c = key[1]
    value = entry_dict[key]
    fout_mat.write(str(row_c) + " " + str(col_c) + " " + str(value) + "\n")

fout_mat.close()