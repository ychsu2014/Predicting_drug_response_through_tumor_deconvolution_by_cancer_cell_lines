#in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/processed_data/oncokb_biomarker_drug_associations_filteredAndExpanded.txt"
#out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/processed_data/oncokb_biomarker_drug_associations_included_final.txt"
#out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/oncoKB/processed_data/oncokb_biomarker_drug_associations_excluded_final.txt"

in_file = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck.txt"
out_file1 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck_included_final.txt"
out_file2 = "/home/CBBI/hsuy1/projects/drugResponse/data/ICIBM_2023/oncoKB/processed_data/oncokb_biomarker_drug_associations_autoCheck_excluded_final.txt"

f = open(in_file)
fout1 = open(out_file1, "w")
fout2 = open(out_file2, "w")

# manually checked that all entries with excluding amino acid all contain "(" and ")"
lines = f.readlines()
header = lines[0]
fout1.write(header)
fout2.write(header.strip("\n") + "\tinclude_AA_change\n")
lines = lines[1:]
for line in lines:
    cols = line.strip("\n").split("\t")
    alt_type = cols[2]
    if "(" in alt_type:
        excl = alt_type.split(" (")[1]
        incl = alt_type.split(" (")[0]
        excl = excl.replace("excluding ", "").replace(")", "").replace(" and ", ", ")
        excls = excl.split(", ")
        for aexc in excls:
            # for excluding these vars in the future, use: gene, incl, cancer type, drug as indexes
            fout2.write("\t".join([cols[0], cols[1], aexc, cols[3], cols[4], incl]) + "\n")
        fout1.write("\t".join([cols[0], cols[1], incl, cols[3], cols[4]]) + "\n")
    else:
        fout1.write(line)

f.close()
fout1.close()
fout2.close()
