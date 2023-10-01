in_ori = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/TCGA_all_mutations_simplifiedColumns_forConvertToGRCh38_removeFailed_v2.txt"
in_convert = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/files_for_hg19_to_GRCh38/hglft_genome_17c0d_9dfb40.bed"
#in_chrom_mapping = "/home/CBBI/hsuy1/projects/drugResponse/data/TCGA/processed_data/mutation/GCA_009914755.4_assembly_report_chromGenBankAccn.txt"
in_TCGA = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns.txt"
out_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38.txt"
out_drop = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns_hg19ConvertedToGRCh38_nomapping.txt"

# accn to chrom
#f_chrom = open(in_chrom_mapping)
#lines_chr = f_chrom.readlines()
#lines_chr = lines_chr[28:]
#accn_to_chrom_dict = {}
#for line in lines_chr:
#    cols = line.strip("\n").split("\t")
#    accn = cols[4]
#    chrom = cols[2]
#    accn_to_chrom_dict[accn] = chrom
#f_chrom.close()

f_ori = open(in_ori)
f_ct = open(in_convert)
lines1 = f_ori.readlines()
lines2 = f_ct.readlines()
print(len(lines1))
print(len(lines2))
ori_to_ct_dict = {}
for i in range(len(lines1)):
    ori_line = lines1[i]
    ct_line = lines2[i]
#    accn = ct_line.split(":")[0]
#    chrom = accn_to_chrom_dict[accn]
#    ct_line = "chr" + chrom + ":" + ct_line.split(":")[1] + "\n"
    ori_to_ct_dict[ori_line] = ct_line
f_ori.close()
f_ct.close()

f = open(in_TCGA)
fout = open(out_file, "w")
fdrop = open(out_drop, "w")
lines = f.readlines()
header = lines[0]
fout.write(header.strip("\n") + "\tGRCh38_chr\tGRCh38_start\tGRCh38_end\n")
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    qline = "chr" + cols[1] + ":" + cols[2] + "-" + cols[3] + "\n"
    if qline in ori_to_ct_dict.keys():
        ct_line = ori_to_ct_dict[qline].strip("\n").replace("chr", "")
        chrom = ct_line.split(":")[0]
        start_pos = ct_line.split(":")[1].split("-")[0]
        end_pos = ct_line.split(":")[1].split("-")[1]
        fout.write(line.strip("\n") + "\t" + "\t".join([chrom, start_pos, end_pos]) + "\n")
    else:
        fdrop.write(line)
f.close()
fout.close()
fdrop.close()
