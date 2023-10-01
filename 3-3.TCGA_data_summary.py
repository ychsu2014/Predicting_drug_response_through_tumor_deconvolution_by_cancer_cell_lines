# cancer types vs. patient num (with all RNA, DNA, and phenotype data)
# the corresponding cell lines among scRNA and COSMIC dataset
in_pheno = "/home/yuching/projects/drugResponse/data/TCGA/raw_data/UCSC_xenabrowser/TCGA_phenotype_denseDataOnlyDownload.tsv"
in_RNA = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/RNA/RNA_total/all_RNA_idx_sID.txt"
in_DNA = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/mutation/TCGA_all_mutations_simplifiedColumns.txt"
out_file = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/TCGA_pheno_RNA_DNA_overlapped_sID_cancerType_list.txt"
out_num = "/home/yuching/projects/drugResponse/data/TCGA/processed_data/TCGA_pheno_RNA_DNA_overlapped_cancerType_pnum.txt"

# sID vs. cancer types
f = open(in_pheno)
lines = f.readlines()
pheno_header = lines[0]
lines = lines[1:]
sID_cancertype_dict = {}
cancertype_sID_dict = {}
for line in lines:
    cols = line.strip("\n").split("\t")
    sID = cols[0]
    cancer_type = cols[3]
    sID_cancertype_dict[sID] = line
    if cancer_type in cancertype_sID_dict.keys():
        old_list = cancertype_sID_dict[cancer_type]
        old_list.append(sID)
        old_list = list(set(old_list))
        cancertype_sID_dict[cancer_type] = old_list
    else:
        cancertype_sID_dict[cancer_type] = [sID]
f.close()
pheno_sID_list = list(set(sID_cancertype_dict.keys()))
print("phenotype data patient num: " + str(len(pheno_sID_list)))

# sID list of RNA dataset
f = open(in_RNA)
lines = f.readlines()
RNA_sID_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    RNA_sID_list.append(cols[1])
f.close()
RNA_sID_list = list(set(RNA_sID_list))
print("RNA data patient num: " + str(len(RNA_sID_list)))

# sID list of DNA dataset
f = open(in_DNA)
lines = f.readlines()
DNA_header = lines[0]
lines = lines[1:]
DNA_sID_list = []
for line in lines:
    cols = line.strip("\n").split("\t")
    DNA_sID_list.append(cols[0])
f.close()
DNA_sID_list = list(set(DNA_sID_list))
print("DNA data patient num: " + str(len(DNA_sID_list)))

fout = open(out_file, "w")
fout.write(pheno_header)
overlapped_pheno_RNA = set(pheno_sID_list).intersection(set(RNA_sID_list))
overlapped_pheno_RNA_DNA = list(overlapped_pheno_RNA.intersection(set(DNA_sID_list)))
for sID in overlapped_pheno_RNA_DNA:
    wline = sID_cancertype_dict[sID]
    fout.write(wline)
fout.close()

fout = open(out_num, "w")
fout.write("cancer_type\tpatient_num\n")
for cancer_type in cancertype_sID_dict.keys():
    sID_list = cancertype_sID_dict[cancer_type]
    filtered_sID_list = list(set(sID_list).intersection(set(overlapped_pheno_RNA_DNA)))
    pnum = len(filtered_sID_list)
    fout.write(cancer_type + "\t" + str(pnum) + "\n")
fout.close()