# remove the headers (3 lines): 1. barcodes, 2. cell lines, 3. pool ID
# keep only  one header: barcode
in_file = "/home/yuching/projects/drugResponse/data/scRNA/raw_data_SCP542/UMIcount_data.txt"
out_file = "/home/yuching/projects/drugResponse/data/scRNA/processed_SCP542/UMI_matrix_data.txt"

f = open(in_file)
fout = open(out_file, "w")

lines = f.readlines()
header = lines[0]
header_cols = header.strip("\n").split("\t")
header_cols = header_cols[1:]
fout.write("GENE")
for col in header_cols:
    fout.write("\t" + col)
fout.write("\n")

lines = lines[3:]
for line in lines:
    fout.write(line) 
    
f.close()
fout.close()
    
