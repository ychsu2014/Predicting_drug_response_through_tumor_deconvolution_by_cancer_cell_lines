in_file = "/home/yuching/projects/drugResponse/data/TCGA/raw_data/UCSC_xenabrowser/mc3.v0.2.8.PUBLIC.xena"
out_file = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns.txt"

f = open(in_file)
fout = open(out_file, "w")

lines = f.readlines()
header = lines[0]
fout.write(header.strip("\n") + "\tAA_change\n" )
lines = lines[1:]
# 0, 6-7, 9-11
for line in lines:
    cols = line.strip("\n").split("\t")
    amino_acid_change = cols[8]
    new_aa = amino_acid_change[2:]
    fout.write("\t".join(cols + [new_aa]) + "\n")

f.close()
fout.close()
