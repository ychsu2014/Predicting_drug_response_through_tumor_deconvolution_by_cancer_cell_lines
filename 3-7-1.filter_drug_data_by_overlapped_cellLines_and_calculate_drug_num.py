# cell column: 0, drug column: 1
in_PRISM = "/home/yuching/projects/drugResponse/data/PRISM/primary-screen-replicate-collapsed-logfold-change.csv"
out_path = "/home/yuching/projects/drugResponse/data/processed_drug_datasets/"
in_PRISM_cell = out_path + "PRISM_scRNA_CCLE_cell_lines_list.txt"
out_PRISM_drug = out_path + "primary-screen-replicate-collapsed-logfold-change_filtered.txt"

def filterForPRISM(inDrug, inCell, outFile, inCellCol1, inCellCol2):
    f1 = open(inDrug)
    f2 = open(inCell)
    # cell list for filtering
    lines2 = f2.readlines()
    cell_list = []
    for line in lines2:
        cols = line.strip("\n").split("\t")
        cell_name = cols[inCellCol2]
        cell_list.append(cell_name)
    cell_list = list(set(cell_list))
    # filter and write to file
    fout = open(outFile, "w")
    lines = f1.readlines()
    header = lines[0]
    fout.write(header)
    drug_list = header.strip("\n").split(",")
    if "" in drug_list:
        drug_list.remove("")
    drug_list = list(set(drug_list))
    lines = lines[1:]
    for line in lines:
        cols = line.strip("\n").split(",")
        cell_name = cols[inCellCol1]
        if cell_name in cell_list:
            fout.write(line)
    print(inDrug)
    print("Number of drugs:")
    print(len(drug_list))
    f1.close()
    f2.close()
    fout.close()
    return(drug_list)

PRISM_drug_list = filterForPRISM(in_PRISM, in_PRISM_cell, out_PRISM_drug, 0, 4)