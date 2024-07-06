import argparse
import os
import re

# author: Martin Kahabka

# structure - can be removed later
# save ID of each analysed patient (c)
# get condition of each patient with saved ID (c)
# counter number of severe/not severe patients (c)
# counter number of  (c)
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

def add_conditions(dict_id, filename):
    num_severe = 0
    num_not_severe = 0
    with open(info_file_path, 'r') as info_file:
        for line in info_file:
            content = line.split("\t")
            # check for header line
            lab_id = content[2]
            condition = content[15]
            # check if lab id correlates to patient
            if lab_id != "Lab ID" and lab_id in id_and_condition:
                dict_id[lab_id] = condition
                # check type of condition
                if condition == "COVID severe":
                    num_severe += 1
                else:
                    num_not_severe += 1
    return (dict_id, num_severe, num_not_severe)

parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-i", "--input_dir", help="Path to the input directory of the filtered vcfs")
parser.add_argument("-p", "--patient_info", help="File with general information of patients")
args = parser.parse_args()

### informationAndData/output_promoter/
# get parameters
input_path = args.input_dir
info_file_path = args.patient_info

# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS FOR PROMTER VCF FILES ---")
counter_pat = 0
id_and_condition = {}
for file in os.listdir(input_path):
    counter_pat += 1
    pattern = r"FO\d*x\d*"
    patient_ID = extract_id(file, pattern)
    id_and_condition[patient_ID] = ""

print("--- SUCCESSFUL: number of files: " + str(counter_pat))
## 
# read conditions from Q001H_sample_preparations_20230803115337.tsv
# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")
id_and_condition, counter_severe, counter_not_severe = add_conditions(id_and_condition, info_file_path)
print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")
