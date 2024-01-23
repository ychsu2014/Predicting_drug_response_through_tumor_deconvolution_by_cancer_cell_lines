# Used the website conversion tool to convert from hg19 to GRCh38 twice.
#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38.txt"
#in_failed = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/hglft_genome_3d7ba_9da170.err.txt"
#out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38_removeFailed.txt"

in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38_removeFailed.txt"
in_failed = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/hglft_genome_116ae_9dd240.err.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38_removeFailed_v2.txt"

f = open(in_failed)
lines = f.readlines()
remove_line_list = []
for line in lines:
    if "#" not in line:
        remove_line_list.append(line)
f.close()

f = open(in_file)
fout = open(out_file, "w")
lines = f.readlines()
for line in lines:
    if line not in remove_line_list:
        fout.write(line)
f.close()
fout.close()
