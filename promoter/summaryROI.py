# author martin kahabka
# use: reads in information from patients and amount of variants per promoter.
#      Then for each promoter collects the amount of variant to either the array of severe covid-19 patients
#      or not severe covid-19 patients
import argparse
import os
import utils.summary_utils as s_utils
    
class ROI_region:
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
            
    def regionToString(self):
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

counter_pat, id_and_condition = s_utils.get_patient_id(input_path, pattern)

print("--- SUCCESSFUL: number of files: " + str(counter_pat))

# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = s_utils.add_conditions(id_and_condition, info_file_path, id_and_condition)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")

regions = {}
# read in all files in input folder and counts variants per ROI
for file_name in os.listdir(input_path):
    print("Working in file: " + file_name)
    # extract patients condition based on file name
    lab_id = s_utils.extract_id(file_name, pattern)
    severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
    
    with open(input_path + "/" + file_name, "r") as roi_file:
        for r in roi_file:
            # read in current promoter
            content = r.split("\t")
            chrom = content[0]
            pos = content[1]
            name = content[2]
            num_variants = content[3].strip() # remove newline char at end
            
            # key for dict
            key = (chrom, pos, name)
            
            # gets promoter by key or adds new promoter to set
            if key in regions:
                current_promoter = regions.get(key)
            else:
                current_promoter = ROI_region(chrom, pos, name)
                regions[key] = current_promoter
            # if patient has severe condition, add to severe array else not severe
            current_promoter.addSum(severe, num_variants)

# save to output file
with open(path_to_output_file, "w") as output_file:
    for current_promoter in regions.values():
        output_file.write(current_promoter.promoterToString() + "\n")
                