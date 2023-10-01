CCLE_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/CCLE/mutation/"
TCGA_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/"
total_out = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/mutation_comparison/TCGA_vs_CCLE/TCGA_CCLE_mut_num.txt"

TCGA_gene_colnum = 6
TCGA_ID_colnum = 0
TCGA_aa_change_colnum = 8
CCLE_gene_colnum = 0
CCLE_ID_colnum = 15
CCLE_aa_change_colnum = 18

import os

CCLE_file_list = os.listdir(CCLE_path)
TCGA_file_list = os.listdir(TCGA_path)
ftotalout = open(total_out, "w")
ftotalout.write("cancer_type\tTCGA_all_mut_num\tCCLE_all_mut_num\tTCGA_not_in_CCLE\tTCGA_in_CCLE\n")

def getTCGAUniMutList(inFile, TCGA_gene_colnum = 6, TCGA_aa_change_colnum = 8, wfile = ftotalout):
    f = open(inFile)
    lines = f.readlines()
    TCGA_header = lines[0]
    lines = lines[1:]   
    #TCGA_uni_mut_list = []
    TCGA_uniMut_line_dict = {}
    for line in lines:
        cols = line.strip("\n").split("\t")
        gene = cols[TCGA_gene_colnum]
        # remove empty entry in hgvsp column
        hgvsp = cols[TCGA_aa_change_colnum]
        if hgvsp != "":
            #TCGA_uni_mut_list.append((gene, hgvsp))
            TCGA_uniMut_line_dict[(gene, hgvsp)] = line
    #TCGA_uni_mut_list = list(set(TCGA_uni_mut_list))
    print("Number of TCGA mutations: ")
    #print(len(TCGA_uni_mut_list))
    print(len(TCGA_uniMut_line_dict.keys()))
    #wfile.write("\t" + str(len(TCGA_uni_mut_list)))
    wfile.write("\t" + str(len(TCGA_uniMut_line_dict.keys())))
    f.close()
    #return(TCGA_uni_mut_list)
    return([TCGA_uniMut_line_dict, TCGA_header])

def CCLEMutCellsDict(inFile, CCLE_gene_colnum = 0, CCLE_aa_change_colnum = 18, CCLE_ID_colnum = 15, wfile = ftotalout):
    f = open(inFile)
    lines = f.readlines()
    lines = lines[1:]
    CCLEmut_cell_dict = {}
    for line in lines:
        cols = line.strip("\n").split(",")
        gene = cols[CCLE_gene_colnum]
        hgvsp = cols[CCLE_aa_change_colnum]
        cell_ID = cols[CCLE_ID_colnum]
        if (gene, hgvsp) in CCLEmut_cell_dict.keys():
            old_list = CCLEmut_cell_dict[(gene, hgvsp)]
            old_list.append(cell_ID)
            old_list = list(set(old_list))
            CCLEmut_cell_dict[(gene, hgvsp)] = old_list
        else:
            CCLEmut_cell_dict[(gene, hgvsp)] = [cell_ID]
    f.close()
    print("CCLE unique mutation: ")
    print(len(CCLEmut_cell_dict.keys()))
    wfile.write("\t" + str(len(CCLEmut_cell_dict.keys())))
    return(CCLEmut_cell_dict)

for CCLE_file in CCLE_file_list:
    in_CCLE = CCLE_path + CCLE_file
    cancer_type = CCLE_file.replace("_filtered_194.txt", "")
    ftotalout.write(cancer_type)
    print(in_CCLE)
    TCGA_file = CCLE_file.replace("_filtered_194.txt", "_GRCh38_7781filtered.txt")
    if TCGA_file in TCGA_file_list:
        in_TCGA = TCGA_path + TCGA_file
        print(in_TCGA)
        out_file = out_path + CCLE_file.replace("_filtered_194.txt", "_TCGA_CCLE_compare.txt")
        out_file2 = out_file.replace(".txt", "_noMatch.txt")
        fout = open(out_file, "w")
        fout2 = open(out_file2, "w")
        #fout.write("TCGA_gene\tTCGA_aa_change\tcell_lines\n")
        #TCGA_uni_mut_list = getTCGAUniMutList(in_TCGA)
        [TCGA_uniMut_line_dict, TCGA_header] = getTCGAUniMutList(in_TCGA)
        fout.write(TCGA_header.strip("\n") + "\tcell_lines\n")
        fout2.write(TCGA_header)
        CCLEmut_cell_dict = CCLEMutCellsDict(in_CCLE)
        TCGA_in_CCLE_mut_num = 0
        TCGA_not_in_CCLE_mut_num = 0
        #for mut in TCGA_uni_mut_list:
        for mut in TCGA_uniMut_line_dict.keys():
            wline = TCGA_uniMut_line_dict[mut]
            if mut in CCLEmut_cell_dict.keys():
                TCGA_in_CCLE_mut_num += 1
                cell_list = CCLEmut_cell_dict[mut]
                fout.write(wline.strip("\n") + "\t" + "\t".join(cell_list) + "\n")
            else:
                TCGA_not_in_CCLE_mut_num += 1
                #fout.write("\n")
                fout2.write(wline)
        print("Number of TCGA mutations not included in CCLE: ")
        print(TCGA_not_in_CCLE_mut_num)
        ftotalout.write("\t" + str(TCGA_not_in_CCLE_mut_num))
        print("Number of TCGA mutations included in CCLE: ")
        print(TCGA_in_CCLE_mut_num)
        ftotalout.write("\t" + str(TCGA_in_CCLE_mut_num) + "\n")
        fout.close()
    else:
        ftotalout.write("\n")
ftotalout.close()
