#in_phe = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7861.txt"
#in_data = "/home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7861.txt"
#out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/"
#out_suffix = "_GRCh38_7861filtered.txt"

in_phe = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/TCGA_filtered_overlapped_sID_cancertype_mapping_7781.txt"
in_data = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_filtered_7781.txt"
out_path = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/"
out_suffix = "_GRCh38_7781filtered.txt"

f = open(in_phe)
lines = f.readlines()
lines = lines[1:]
cancertype_sID_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    sID = cols[0]
    sc_cancer_type = cols[2]
    if sc_cancer_type in cancertype_sID_dict.keys():
        old_list = cancertype_sID_dict[sc_cancer_type]
        old_list.append(sID)
        old_list = list(set(old_list))
        cancertype_sID_dict[sc_cancer_type] = old_list
    else:
        cancertype_sID_dict[sc_cancer_type] = [sID]
f.close()

f = open(in_data)
lines = f.readlines()
header = lines[0]
lines = lines[1:]
fout1 = open(out_file, "w")
fout1.write(header.strip("\n") + "\tscRNA_cancer_type\n")
for cancer_type in cancertype_sID_dict.keys():
    sID_list = cancertype_sID_dict[cancer_type]
    w_file_prefix = cancer_type
    if " " in w_file_prefix:
        w_file_prefix = w_file_prefix.replace(" ", "")
    if "/" in w_file_prefix:
        w_file_prefix = w_file_prefix.replace("/", "")
    fout = open(out_path + w_file_prefix + out_suffix, "w")
    fout.write(header)
    for line in lines:
        cols = line.strip("\n").split("\t")
        sID = cols[0]
        if sID in sID_list:
            fout.write(line)
            fout1.write(line.strip("\n") + "\t" + cancer_type + "\n")
    fout.close()
fout1.close()
f.close()
