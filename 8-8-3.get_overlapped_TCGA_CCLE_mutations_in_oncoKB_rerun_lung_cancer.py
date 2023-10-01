# only include HGVSP and oncogenic mutations
in_keep = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck_included_final.txt"
in_drop = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck_excluded_final.txt"
#in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/"

in_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/lung_cancer_all_oncoKB_19/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE_in_oncoKB/lung_cancer_all_oncoKB_19/"

import re
import os

in_files_temp = os.listdir(in_path)
in_files = []
for in_file in in_files_temp:
    if "_TCGA_CCLE_compare.txt" in in_file:
        in_files.append(in_file)
in_files = list(set(in_files))

pattern_full = "[A-Z]{1}\d+[A-Z]{1}"
pattern_all_aa = "[A-Z]{1}\d+"

# parse oncoKB data
def parseMutToDicts(inFile):
    pattern_all_aa = "[A-Z]{1}\d+"
    f = open(inFile)
    lines = f.readlines()
    lines = lines[1:]
    mutLineDict = {}
    allAASubList = []
    oncoMutList = []
    wildtypeList = []
    otherList = []
    for line in lines:
        cols = line.strip("\n").split("\t")
        gene = cols[1]
        alt = cols[2]
        # alt: aa changes (hgvsp), Oncogenic Mutations, Wildtype
        if (gene, alt) in mutLineDict.keys():
            oldList = mutLineDict[(gene, alt)]
            oldList.append(line)
            mutLineDict[(gene, alt)] = oldList
        else:
            mutLineDict[(gene, alt)] = [line]
        # special handle 1: oncogenic mutations
        if alt.lower() == "oncogenic mutations":
            oncoMutList.append(gene)
            oncoMutList = list(set(oncoMutList))
        # special handle 2: wildtype
        elif alt.lower() == "wildtype":
            wildtypeList.append(gene)
            wildtypeList = list(set(wildtypeList))
        elif re.match(pattern_all_aa, alt) != None:
            # special handle 3: all aa substitution
            if re.match(pattern_all_aa, alt)[0] == alt:
                allAASubList.append([gene, alt])
            else:
                otherList.append([gene, alt])
        else:
            otherList.append([gene, alt])
    f.close()
    return([allAASubList, oncoMutList, wildtypeList, otherList, mutLineDict])
            
# mutations to exclude
# onco_drop and wild_drop are empty (did not do anything to handle)
[all_aa_drop, onco_drop, wild_drop, other_drop, drop_mut_line_dict] = parseMutToDicts(in_drop)
#print(all_aa_drop)
#print(onco_drop)
#print(wild_drop)
#print(other_drop)

# mutations to keep
[all_aa_keep, onco_keep, wild_keep, other_keep, keep_mut_line_dict] = parseMutToDicts(in_keep)
#print(all_aa_keep)
#print(onco_keep)
#print(wild_keep)
#print(other_keep)

def getLineAndWrite(inGeneAlt, inDict, inTCGALine, inFile, checkAllAALines):
    allAAIdx = False
    gotLineIdx = False
    if checkAllAALines == True:
        #print("check AA")
        pattern2 = "[A-Z]{1}\d+"
        alt = inGeneAlt[1]
        mResult = re.match(pattern2, alt)
        if mResult != None:
            alt = mResult[0]
            if (inGeneAlt[0], alt) in inDict.keys():
                wlines = inDict[(inGeneAlt[0], alt)]
                for wline in wlines:
                    inFile.write(inTCGALine.strip("\n") + "\t" + wline)
                allAAIdx = True
    if inGeneAlt in inDict.keys():
        #print("checked all other vars.")
        wlines = inDict[inGeneAlt]
        for wline in wlines:
            inFile.write(inTCGALine.strip("\n") + "\t" + wline)
        gotLineIdx = True
    if allAAIdx == False and gotLineIdx == False:
        print("Please check")
        wline = ""
        inFile.write(inTCGALine.strip("\n") + "\t" + wline + "\n")

#def wildtypeGetLineAndWrite(inGeneAlt, inDict, inFile, inSID):
#    if inGeneAlt in inDict.keys():
#        wlines = inDict[inGeneAlt]
#        for wline in wlines:
#            inFile.write(inSID + "\t" + wline)

def ifInAASubLists(inGene, inAlt, inOtherList, inAllAASubList, matchAllAASubIdx):
    outInListIdx = False
    if matchAllAASubIdx == True:
        pattern2 = "[A-Z]{1}\d+"
        mResult = re.match(pattern2, inAlt)
        if mResult != None:
            if [inGene, mResult[0]] in inAllAASubList:
                outInListIdx = True
    if [inGene, inAlt] in inOtherList:
        outInListIdx = True
    return(outInListIdx)

#def writeWildtype(inFile, outWildType, inWildKeepList = wild_keep, inKeepMutLineDict = keep_mut_line_dict):
#    f_in = open(inFile)
#    lines = f_in.readlines()
#    header = lines[0]
#    lines = lines[1:]
#    # this dict is for checking wildtype
#    sID_gene_dict = {}
#    fwildtype = open(outWildType, "w")
#    for line in lines:
#        cols = line.strip("\n").split("\t")
#        sID = cols[0]
#        gene = cols[6]
#        if sID in sID_gene_dict.keys():
#            old_list = sID_gene_dict[sID]
#            old_list.append(gene)
#            old_list = list(set(old_list))
#            sID_gene_dict[sID] = old_list
#        else:
#            sID_gene_dict[sID] = [gene]
#    # write to file
#    for sID in sID_gene_dict.keys():
#        gene_list = sID_gene_dict[sID]
#        wildtype_genes = list(set(inWildKeepList).difference(set(gene_list)))
#        for gene in wildtype_genes:
#            # getLineAndWrite((gene, "Wildtype"), inKeepMutLineDict, "", fwildtype, False)
#            wildtypeGetLineAndWrite((gene, "Wildtype"), inKeepMutLineDict, fwildtype, sID)
#    fwildtype.close()

def writeToFile(inFile, outFile, outDrop, outExclude, inOtherKeep = other_keep, inAllaaKeep = all_aa_keep, inWildKeep = wild_keep, inKeepMutLineDict = keep_mut_line_dict, inOtherDrop = other_drop, inAllaaDrop = all_aa_drop, inDropMutLineDict = drop_mut_line_dict, inOncoKeep = onco_keep):
    pattern_full = "[A-Z]{1}\d+[A-Z]{1}"
    f_in = open(inFile)
    lines = f_in.readlines()
    header = lines[0]
    lines = lines[1:]
    fout = open(outFile, "w")
    fdrop = open(outDrop, "w")
    fexclude = open(outExclude, "w")
    for line in lines:
        cols = line.strip("\n").split("\t")
        gene = cols[6]
        alt = cols[12]
        if "splice" in alt or "ext" in alt or "delins" in alt or "del" in alt or "dup" in alt or "ins" in alt or "fs" in alt:
            # exclude mutations
            exclude_idx = ifInAASubLists(gene, alt, inOtherDrop, inAllaaDrop, False)
            if exclude_idx == True:
                getLineAndWrite((gene, alt), inDropMutLineDict, line, fexclude, False)
            # include mutations
            include_idx = ifInAASubLists(gene, alt, inOtherKeep, inAllaaKeep, False)
            if include_idx == True:
                getLineAndWrite((gene, alt), inKeepMutLineDict, line, fout, False)
            if gene in inOncoKeep:
                getLineAndWrite((gene, "Oncogenic Mutations"), inKeepMutLineDict, line, fout, False)
                onco_keep_idx = True
            else:
                onco_keep_idx = False
            # remove the line that was included in the wildtype output file
            if gene in inWildKeep:
                #print(gene)
                #print(alt)
                getLineAndWrite((gene, "Wildtype"), inKeepMutLineDict, line, fexclude, False)
                wild_keep_idx = True
            else:
                wild_keep_idx = False
            if exclude_idx == False and include_idx == False and onco_keep_idx == False and wild_keep_idx == False:
                fdrop.write(line)
        elif "?" in alt:
            fdrop.write(line)
        elif re.match(pattern_full, alt) != None:
            # exclude mutations
            exclude_idx = ifInAASubLists(gene, alt, inOtherDrop, inAllaaDrop, True)
            if exclude_idx == True:
                getLineAndWrite((gene, alt), inDropMutLineDict, line, fexclude, True)
            # some mutations may be excluded and also included for different drugs
            include_idx = ifInAASubLists(gene, alt, inOtherKeep, inAllaaKeep, True)
            if include_idx == True:
                getLineAndWrite((gene, alt), inKeepMutLineDict, line, fout, True)
            if gene in inOncoKeep:
                getLineAndWrite((gene, "Oncogenic Mutations"), inKeepMutLineDict, line, fout, False)
                onco_keep_idx = True
            else:
                onco_keep_idx = False
            # remove the line that was included in the wildtype output file
            if gene in inWildKeep:
                getLineAndWrite((gene, "Wildtype"), inKeepMutLineDict, line, fexclude, False)
                wild_keep_idx = True
            else:
                wild_keep_idx = False
            if exclude_idx == False and include_idx == False and onco_keep_idx == False and wild_keep_idx == False:
                fdrop.write(line)
        else:
             fdrop.write(line)
    f_in.close()
    fout.close()
    fdrop.close()
    fexclude.close()

for in_file in in_files:
    out_file = out_path + in_file.replace(".txt", "_in_oncoKB.txt")
    out_drop = out_path + in_file.replace(".txt", "_oncoKB_dropped.txt")
    out_exclude = out_path + in_file.replace(".txt", "_oncoKB_excluded.txt")
    out_wildtype = out_path + in_file.replace(".txt", "_oncoKB_wildtype.txt")
    in_file = in_path + in_file
    #writeWildtype(in_file, out_wildtype)
    writeToFile(in_file, out_file, out_drop, out_exclude)
