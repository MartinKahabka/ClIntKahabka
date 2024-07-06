import argparse
import os
import re

# author: Martin Kahabka

# this file should input the output if "runVariantProm.sh" and return a tsv as described on git
# necessary files
# 1. from output files of "runVariantProm.sh" get
#   1.1 chr/pos of variants
#   1.2 number of patients with variant
#   1.3 ID of patients
# 3. from Q001H_sample_preparations_20230803115337.tsv get
#   3.1 condition (severe/not severe) of patients

# structure
# save ID of each analysed patient
# get condition of each patient with saved ID
# counter number of severe/not severe patients
# counter number of 
# for each output file get
#   ID of patient -> condition of patient
#   for each variant do
#       is variant already saved? -> create new dict entry
#       increase correlating counter
#       save ID to correlation array
# iterate over variants and fill number of patients where variant was not found (see number severe/not severe patients)
# save file in described format (see git)

# gets a string and a regex pattern that describes the unique ID of the patient
def extract_id(name, pattern) -> str:
    id = re.search(pattern, name).group()
    return id


parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-i", "--input_dir", help="Path to the input directory")
args = parser.parse_args()

### informationAndData/output_promoter/
# get parameters
input_path = args.input_dir

# get IDs of patients
ids = []
for file in os.listdir(input_path):
    pattern = r"FO\d*x\d*"
    patient_ID = extract_id(file, pattern)
    ids.append(patient_ID)

print(ids)