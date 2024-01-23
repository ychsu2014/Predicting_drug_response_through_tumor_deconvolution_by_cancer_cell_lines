in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38.txt"

f = open(in_file)
fout = open(out_file, "w")
lines = f.readlines()
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    chrom = cols[1]
    start_pos = cols[2]
    end_pos = cols[3]
    fout.write("chr" + chrom + ":" + start_pos + "-" + end_pos + "\n")
f.close()
fout.close()
    
