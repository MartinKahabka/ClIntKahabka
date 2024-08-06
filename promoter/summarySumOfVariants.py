import argparse
import re
import os
# for each promoter: collect number of variants for severe and not severe
    # get condition of patient
    # get sum of promoter variants
    # add sums to dict
    
class promoter:
    def __init__(self, chrom, pos, name):
        self.chr = chrom
        self.pos = pos
        self.name = name
        self.severeCov = []
        self.notSevereCov = []
        
    def addSum(self, severe, num):
        if severe:
            self.severeCov.append(num)
        else:
            self.notSevereCov.append(num)
            
    def promoterToString(self):
        header = "\t".join([self.chr, self.pos, self.name])
        pos = "\t".join(self.severeCov)
        neg = "\t".join(self.notSevereCov)
        print(self.severeCov)
        print(repr(pos))
        return "\n".join([header, pos, neg])
        
# gets a string and a regex pattern that describes the unique ID of the patient
def extract_id(name, pattern) -> str:
    try:
        id = re.search(pattern, name).group()   
    except AttributeError:
        id = "nullID"
    return id

def add_conditions(dict_id, filename):
    num_severe = 0
    num_not_severe = 0
    with open(filename, 'r') as info_file:
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

            
print("--- START PROGRAM SUMMARYSUMOFVARIANTS.PY ---")

parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-i", "--input_dir", help="Path to the input directory of the sum of variants per promoter per patient")
parser.add_argument("-p", "--patient_info", help="File with general information of patients")
parser.add_argument("-o", "--output_dir", help="Path of output file")
args = parser.parse_args()

### informationAndData/output_promoter/
# get parameters
input_path = args.input_dir
info_file_path = args.patient_info
output_file_path = args.output_dir

pattern = r"FO\d*x\d*"
name_output_file = "summary_sum_variants.tsv"
path_to_output_file = os.path.join(output_file_path, name_output_file)

# print args
print(input_path)
print(info_file_path)
print(output_file_path)


# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS FOR PROMTER VCF FILES ---")

counter_pat = 0
id_and_condition = {}
for file in os.listdir(input_path):
    counter_pat += 1
    patient_ID = extract_id(file, pattern)
    if patient_ID != "nullID":
        id_and_condition[patient_ID] = ""

print("--- SUCCESSFUL: number of files: " + str(counter_pat))

# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = add_conditions(id_and_condition, info_file_path)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")

sum_promoter = {}
for file_name in os.listdir(input_path):
    print("Working in file: " + file_name)
    # get lab id
    lab_id = extract_id(file_name, pattern)
    severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
    
    with open(input_path + "/" + file_name, "r") as promoter_file:
        for p in promoter_file:
            # read in promoter
            content = p.split("\t")
            chrom = content[0]
            pos = content[1]
            name = content[2]
            num_variants = content[3].strip() # remove newline char at end
            
            # key for dict
            key = (chrom, pos, name)
            
            if key in sum_promoter:
                # if patient has severe condition
                current_promoter = sum_promoter.get(key)
            else:
                current_promoter = promoter(chrom, pos, name)
                sum_promoter[key] = current_promoter
            # add to promoter
            current_promoter.addSum(severe, num_variants)

# save to output file
with open(path_to_output_file, "w") as output_file:
    for current_promoter in sum_promoter.values():
        output_file.write(current_promoter.promoterToString() + "\n")
                