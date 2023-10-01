in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/raw_data/oncokb_biomarker_drug_associations.tsv"
# keep only point mutations and reformat one mutation one drug per line (originally, there are multiple mutations one drug per line)
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_manualCheck.txt"
out_file3 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_drop.txt"

import re

f = open(in_file)
fout1 = open(out_file1, "w")
fout2 = open(out_file2, "w")
fout3 = open(out_file3, "w")

lines = f.readlines()
header = lines[0]
fout1.write(header)
fout2.write(header)
fout3.write(header)
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    alt_types = cols[2]
    if "(" not in  alt_types and ")" not in alt_types:
        alt_types = alt_types.split(", ")
    else:
        # manually checked that all entries with "(" and ")" contain only one type of variants
        alt_types = [alt_types]
    for alt_type in alt_types:
        out_line = "\t".join([cols[0], cols[1], alt_type, cols[3], cols[4]]) + "\n"
        alt_type_lower = alt_type.lower()
        if "fusion" in alt_type_lower or "deletion" in alt_type_lower or "duplication" in alt_type_lower or "insertions" in alt_type_lower or "amplification" in alt_type_lower or alt_type_lower == "Internal tandem duplication".lower():
            fout2.write(out_line)
        elif alt_type_lower == "Tumor Mutational Burden-High".lower() or alt_type_lower == "Microsatellite Instability-High".lower():
            fout3.write(out_line)
        elif (re.match(".*del", alt_type_lower) != None) or (re.match(".*dup", alt_type_lower) != None) or (re.match(".*ins.*", alt_type_lower) != None) or (re.match(".*delins.*", alt_type_lower) != None) or (re.match(".*splice", alt_type_lower) != None):
            fout1.write(out_line)
        elif ("Oncogenic Mutations".lower() not in alt_type_lower) and (re.match(".*mutations.*", alt_type_lower) != None):
            fout2.write(out_line)
        else:
            fout1.write(out_line)

f.close()
fout1.close()
fout2.close()
fout3.close()
