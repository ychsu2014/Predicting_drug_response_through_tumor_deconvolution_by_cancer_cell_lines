in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/TCGA_CCLE_in_oncoKB_uni_entry_with_duplicated_vars_v2.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/TCGA_CCLE_in_oncoKB_uni_entry_with_duplicated_vars_v2_oneDrugOneLine.txt"

f = open(in_file)
lines = f.readlines()
fout = open(out_file, "w")

# inLines: lines without header; inDrugcol: the number of the column of the oncoKB drugs
# keep the lines with empty drug columns
def oneDrugOneLine(inLines, inDrugcol):
    outLines = []
    for line in inLines:
        cols = line.strip("\n").split("\t")
        drugs = cols[inDrugcol]
        drugList = drugs.split(", ")
        for adrug in drugList:
            outLine = cols[:inDrugcol+1] + [adrug] + cols[inDrugcol+1:]
            outLine = "\t".join(outLine) + "\n"
            outLines.append(outLine)
    return(outLines)

out_lines = oneDrugOneLine(lines, 21)
for out_line in out_lines:
    fout.write(out_line)

f.close()
fout.close()
