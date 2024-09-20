# author martin kahabka
# use: reads in information from patients and amount of variants per promoter.
#      Then for each promoter collects the amount of variant to either the array of severe covid-19 patients
#      or not severe covid-19 patients
import argparse
import re
import os

    
class promoter:
    def __init__(self, chrom, pos, name):
        self.chr = chrom
        self.pos = pos
        self.name = name
        self.severeCov = []
        self.notSevereCov = []
        
    def addSum(self, severe, num):
        """ adds the amount of variants of one patients either to the severe or not severe array

        Args:
            severe (boolean): true if severe condition, false not not severe condition
            num (integer): amount of variants in promoter region of one patient
        """
        if severe:
            self.severeCov.append(num)
        else:
            self.notSevereCov.append(num)
            
    def promoterToString(self):
        """ creates string representation based on the current values of this class

        Returns:
            string: string representation of promoter
        """
        sum_variants_severe = str(sum([int(i) for i in self.severeCov]))
        sum_variants_not = str(sum([int(i) for i in self.notSevereCov]))
        header = "\t".join([self.chr, self.pos, self.name, sum_variants_severe, sum_variants_not])
        pos = "\t".join(self.severeCov)
        neg = "\t".join(self.notSevereCov)
        return "\n".join([header, pos, neg])
        
def extract_id(name, pattern) -> str:
    """ extracts id of the patient from name of file name using a given pattern

    Args:
        name (string): name of patient file
        pattern (_type_): pattern to match string

    Returns:
        str: id of patient, "nullID" if pattern cannot be matched to patient file
    """
    try:
        id = re.search(pattern, name).group()
        return id   
    except AttributeError:
        id = "nullID"
        return id

def get_conditions(dict_id, filename):
    """ finds the patients condition, counts the number of severe and not severe patients

    Args:
        dict_id (dict[str:str]): dict with id of patients
        filename (_type_): path to file with metadata of patients

    Returns:
        (dict[str:str], integer, integer): updated dict with conditions of patients, amount of severe patients/not severe patients
    """
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

def get_patient_id(variant_file_path, pattern_id):
    """ reads in the file names from the input folder and creates a dict with the patient ids. Counts the number of 

    Args:
        variant_file_path (string): path to patients .vcf file
        pattern_id (string): pattern to match

    Returns:
        (integer, dict[string : string]): number of patients and dict with patients id
    """
    counter_patients = 0
    id_and_condition = {}
    for metafile in os.listdir(variant_file_path):
        counter_patients += 1
        patient_ID = extract_id(metafile, pattern_id)
        if patient_ID != "nullID":
            id_and_condition[patient_ID] = ""
    
    return counter_patients, id_and_condition

            
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
path_to_output_file = args.output_dir

pattern = r"FO\d*x\d*"

# print args
print(input_path)
print(info_file_path)
print(path_to_output_file)


# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS FOR PROMTER VCF FILES ---")

counter_pat, id_and_condition = get_patient_id(input_path, pattern)

print("--- SUCCESSFUL: number of files: " + str(counter_pat))

# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = get_conditions(id_and_condition, info_file_path)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")

set_of_promoter = {}
# read in patient data and create arrays with amount of variants per promoter
for file_name in os.listdir(input_path):
    print("Working in file: " + file_name)
    # get lab id and condition of patient
    lab_id = extract_id(file_name, pattern)
    severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
    
    with open(input_path + "/" + file_name, "r") as promoter_file:
        for p in promoter_file:
            # read in current promoter
            content = p.split("\t")
            chrom = content[0]
            pos = content[1]
            name = content[2]
            num_variants = content[3].strip() # remove newline char at end
            
            # key for dict
            key = (chrom, pos, name)
            
            # gets promoter by key or adds new promoter to set
            if key in set_of_promoter:
                current_promoter = set_of_promoter.get(key)
            else:
                current_promoter = promoter(chrom, pos, name)
                set_of_promoter[key] = current_promoter
            # if patient has severe condition, add to severe array else not severe
            current_promoter.addSum(severe, num_variants)

# save to output file
with open(path_to_output_file, "w") as output_file:
    for current_promoter in set_of_promoter.values():
        output_file.write(current_promoter.promoterToString() + "\n")
                