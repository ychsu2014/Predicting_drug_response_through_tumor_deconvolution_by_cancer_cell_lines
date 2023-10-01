in_file = "/home/yuching/projects/drugResponse/data/PRISM/primary-screen-replicate-collapsed-treatment-info.csv"
out_file = "/home/yuching/projects/drugResponse/data/processed_drug_datasets/PRISM_primary_treatment_info.txt"

import pandas as pd

df = pd.read_csv(in_file)
df.to_csv(out_file, sep = "\t", index = False)
