in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/"

import os

all_files = os.listdir(in_path)
inc_files = []
for afile in all_files:
    if "_in_oncoKB.txt" in afile:
        inc_files.append(afile)

def removeExcludedVar(inFile, inExc, outFile, outExc):
    f_exc = open(inExc)
    lines = f_exc.readlines()
    excluded_lines = []
    for line in lines:
        # do not need to handle wildtype
        if "Wildtype" not in line:
            cols = line.strip("\n").split("\t")
            excluded_lines.append(line)
    f_exc.close()
    # write to file
    f_inc = open(inFile)
    fout = open(outFile, "w")
    fdrop = open(outExc, "w")
    lines = f_inc.readlines()
    for line in lines:
        if line in excluded_lines:
            fdrop.write(line)
        else:
            fout.write(line)
    f_inc.close()
    fout.close()
    fdrop.close()

for inc_file in inc_files:
    inc_file = in_path + inc_file
    exc_file = inc_file.replace("_in_oncoKB.txt", "_oncoKB_excluded.txt")
    out_file = inc_file.replace("_in_oncoKB.txt", "_in_oncoKB_final.txt")
    drop_file = inc_file.replace("_in_oncoKB.txt", "_in_oncoKB_excluded.txt")
    removeExcludedVar(inc_file, exc_file, out_file, drop_file)
