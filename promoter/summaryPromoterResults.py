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

# input: dict, str
# dict_id: dict with lab ids as keys
# filename: path to file that contains the condition of the patients
# output: (dict, int, int)
# dict_id: input dict with added conditions
# num_severe: number of patients with condition severe
# num_not_severe: number of patients with other condition
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

# gets dict and fills up the negative counter (num of patients that do not have the mutation) 
# of the dictionary according to the counter of num_severe and num_not_severe
def fill_up_dict(dict_variant, num_severe, num_not_severe) -> dict:
    # fill up posNegative and negNegative
    for variant in dict_variant:
        value = dict_variant.get(variant)
        value[1] = num_severe - value[0]
        value[3] = num_not_severe - value[2]
        dict_variant[variant] = value
    return dict_variant

# gets the information about a variant and returns a string containing the information
def variant_to_string(name, infos) -> str:
    str_infos = ""
    for i in infos:
        str_infos += "\t" + str(i)
    str_name = name[0] + "\t" + str(name[1])
    return str_name + str_infos

def comments_file(filename, input_file_path, num_severe, num_not_severe) -> str:
    # extract unique name from pipeline
    unique_name = os.path.split(input_file_path)[-1]
    # build comment
    comment = ""
    comment += "# unique identifier: " + unique_name + "\n"
    comment += "# input from " + input_file_path + "\n"
    comment += "# output to " + filename + "\n"
    comment += "# number of patients: " + str(num_severe + num_not_severe) + "\n"
    comment += "# number of severe/not severe patients: " + str(num_severe) + "/" + str(num_not_severe) + "\n"
    comment += "# pos/neg corresponds to number of patients where variant is there/not there" + "\n"
    comment += "# severe/NotSevere corresponds to number of patients with have severe/not severe condition" + "\n"
    comment += "chromosome\tposition\tposSevere\tnegSevere\tposNotSevere\tnegNotSevere" + "\n"
    return comment

parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-i", "--input_dir", help="Path to the input directory of the filtered vcfs")
parser.add_argument("-p", "--patient_info", help="File with general information of patients")
parser.add_argument("-o", "--output_dir", help="Path of output file")
args = parser.parse_args()

### informationAndData/output_promoter/
# get parameters
input_path = args.input_dir
info_file_path = args.patient_info
output_file_path = args.output_dir
pattern = r"FO\d*x\d*"

# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS FOR PROMTER VCF FILES ---")

counter_pat = 0
id_and_condition = {}
for file in os.listdir(input_path):
    counter_pat += 1
    patient_ID = extract_id(file, pattern)
    id_and_condition[patient_ID] = ""

print("--- SUCCESSFUL: number of files: " + str(counter_pat))

# read conditions from Q001H_sample_preparations_20230803115337.tsv
# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = add_conditions(id_and_condition, info_file_path)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")


# dict entry: key = chr/pos, value = (severePos, severeNeg, notSeverePos, not SevereNeg)
# define
# mutation [(str)]
print("--- READ IN VARIANTS FROM PROMTER VCF FILES ---")
# read in all variants from each file
variant_information = {}
for file_name in os.listdir(input_path):
    print("Working in file: " + file_name)
    # get lab id
    lab_id = extract_id(file_name, pattern)
    with open(input_path + "/" + file_name, 'r') as file:
        for line in file:
            # check for comments
            if line[0] != '#':
                content = line.split('\t')
                # get relevant values
                chrom = content[0]
                pos = int(content[1])
                severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
                # check if variant is in dict
                if (chrom, pos) in variant_information:
                    # update dict
                    value = variant_information.get((chrom, pos))
                    if severe:
                        value[0] += 1
                    else:
                        value[2] += 1
                else:
                    # create entry
                    if severe:
                        variant_information[(chrom, pos)] = [1, 0, 0, 0]
                    else:
                        variant_information[(chrom, pos)] = [0, 0, 1, 0]
                        
print("--- SUCCESFUL ---")
                    
print("--- SAVE FILE TO " + output_file_path + " ---")
# fill up dict
variant_information = fill_up_dict(variant_information, counter_severe, counter_not_severe)
# save to file
with open(output_file_path, 'w') as output_file:
    c = comments_file(output_file_path, input_path, counter_severe, counter_not_severe)
    output_file.write(c)
    for variant in variant_information:
        s = variant_to_string(variant, variant_information[variant])
        output_file.write(s + "\n")
print("--- SUCCESSFUL -> PROGRAMM TERMINATES ---")