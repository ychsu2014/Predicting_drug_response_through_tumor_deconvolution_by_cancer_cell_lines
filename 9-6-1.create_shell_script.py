out_file = "16-0-1.create_scripts.sh"

cancer_type_list = []
cancer_type_list.append("BileDuctCancer")
cancer_type_list.append("BladderCancer")
cancer_type_list.append("BrainCancer")
cancer_type_list.append("BreastCancer")
cancer_type_list.append("ColonColorectalCancer")
cancer_type_list.append("EndometrialUterineCancer")
cancer_type_list.append("EsophagealCancer")
cancer_type_list.append("GastricCancer")
cancer_type_list.append("HeadandNeckCancer")
cancer_type_list.append("KidneyCancer")
cancer_type_list.append("LiverCancer")
cancer_type_list.append("LungCancer")
cancer_type_list.append("OvarianCancer")
cancer_type_list.append("PancreaticCancer")
cancer_type_list.append("ProstateCancer")
cancer_type_list.append("Sarcoma")
cancer_type_list.append("SkinCancer")
cancer_type_list.append("ThyroidCancer")

partNum_CGCnum_dict = {}
partNum_CGCnum_dict[1] = (0,100)
partNum_CGCnum_dict[2] = (100,200)
partNum_CGCnum_dict[3] = (200,300)
partNum_CGCnum_dict[4] = (300,400)
partNum_CGCnum_dict[5] = (400,500)
partNum_CGCnum_dict[6] = (500,600)
partNum_CGCnum_dict[7] = (600,800)

fout = open(out_file, "w")

fout.write("# 737 CGC genes\n")
fout.write("# cancer_type, part_num, CGC_start_num, CGC_end_num\n")
# python3 16-0_sub.create_scripts.py "BrainCancer" 1 0 100
for cancer_type in cancer_type_list:
	for part_num in partNum_CGCnum_dict.keys():
		(CGC_start_num, CGC_end_num) = partNum_CGCnum_dict[part_num]
		fout.write("python3 16-0-1_sub.create_scripts.py '" + cancer_type + "' " + str(part_num) + " " + str(CGC_start_num) + " " + str(CGC_end_num) + "\n")
fout.close()
