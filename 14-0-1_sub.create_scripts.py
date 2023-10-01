import sys

cancer_type = sys.argv[1]
part_num = sys.argv[2]
CGC_start_num = sys.argv[3]
CGC_end_num = sys.argv[4]

template_file = "16.drugResponse_exp_high_low_comparison_v2.py"
out_file = "16.drugResponse_exp_high_low_comparison_" + cancer_type + "_part" + part_num + "_v2.py"

fout = open(out_file, "w")
fout.write("in_cancer_type = '" + str(cancer_type) + "'\n")
fout.write("CGC_start_num = " + str(CGC_start_num) + "\n")
fout.write("CGC_end_num = " + str(CGC_end_num) + "\n")
fout.write("part_num = " + str(part_num) + "\n")

f = open(template_file)
lines = f.readlines()
for line in lines:
	fout.write(line)

f.close()
fout.close()
