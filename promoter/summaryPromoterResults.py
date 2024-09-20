# author: Martin Kahabka
import argparse
import os
import re

# gets a string and a regex pattern that describes the unique ID of the patient
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
    except AttributeError:
        id = "nullID"
    return id

# input: dict, str
# dict_id: dict with lab ids as keys
# filename: path to file that contains the condition of the patients
# output: (dict, int, int)
# dict_id: input dict with added conditions
# num_severe: number of patients with condition severe
# num_not_severe: number of patients with other condition
def add_conditions(dict_id, filename):
    """ matches conditions (severe/not severe) to the the patients ids

    Args:
        dict_id ( dict[string : string] ): input dict without conditions
        filename (string): path to file that contains the condition of the patients

    Returns:
        ( dict[string : string], integer, integer): input dict with added conditions \n
                                                    number of patients with severe condition \n
                                                    number of patients with not severe condition
    """
    num_severe = 0
    num_not_severe = 0
    # open metadata fole
    with open(filename, 'r') as info_file:
        for line in info_file:
            # get information from line
            content = line.split("\t")
            lab_id = content[2]
            condition = content[15]
            
            # check for header line and check if lab id correlates to patient
            if lab_id != "Lab ID" and lab_id in id_and_condition:
                dict_id[lab_id] = condition
                # add to counter
                if condition == "COVID severe":
                    num_severe += 1
                else:
                    num_not_severe += 1
    return (dict_id, num_severe, num_not_severe)

# gets dict and fills up the negative counter (num of patients that do not have the mutation) 
# of the dictionary according to the counter of num_severe and num_not_severe
# fill up posNegative and negNegative
def fill_up_dict(dict_variant, num_severe, num_not_severe) -> dict:
    """ fills in the number of patients with not variants (posNegative and negNegativ) 
        based on the total number of patients and the number of patients with variants

    Args:
        dict_variant ( dict[string : array[integer]] ): dict with variant information
        num_severe (_type_): number of patients with severe condition
        num_not_severe (_type_): number of patients without severe condition

    Returns:
        dict: updated dict containing the variant information 
    """
    # iterates over each found variant
    for variant in dict_variant:
        value = dict_variant.get(variant)
        # calculates the number of patients without this variant
        value[1] = num_severe - value[0]
        value[3] = num_not_severe - value[2]
        dict_variant[variant] = value
    return dict_variant

# gets the information about a variant and returns a string containing the information
def variant_to_string(name, infos) -> str:
    """ gets key and value of a single variant from dict and returns the information as a string

    Args:
        name ( (string, integer) ): chromosome and position of variant on genome
        infos ( array[integer] ): array containing the amount of variants per category

    Returns:
        str: string representation of the variants information
    """
    # amount of variants per category
    str_infos = ""
    for i in infos:
        str_infos += "\t" + str(i)
    # chrom and pos of variant
    str_name = name[0] + "\t" + str(name[1])
    return str_name + str_infos

def comments_file(filename, input_file_path, num_severe, num_not_severe, numVariants) -> str:
    """ creates comment for file containing useful information

    Args:
        filename (string): name of output file
        input_file_path (string): name of input file
        num_severe (integer): number of patients with severe condition
        num_not_severe (integer): number of patients without severe condition
        numVariants (integer): total number of found variants

    Returns:
        str: _description_
    """
    # extract unique name from pipeline
    unique_name = os.path.split(input_file_path)[-1]
    # build comment
    comment = ""
    comment += "# unique identifier: " + unique_name + "\n"
    comment += "# input from " + input_file_path + "\n"
    comment += "# output to " + filename + "\n"
    comment += "# number of patients: " + str(num_severe + num_not_severe) + "\n"
    comment += "# number of severe/not severe patients: " + str(num_severe) + "/" + str(num_not_severe) + "\n"
    comment += "# number of variants found: " + str(numVariants) + "\n"
    comment += "# pos/neg corresponds to number of patients where variant is there/not there" + "\n"
    comment += "# severe/NotSevere corresponds to number of patients with have severe/not severe condition" + "\n"
    comment += "chromosome\tposition\tposSevere\tnegSevere\tposNotSevere\tnegNotSevere" + "\n"
    return comment

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

def add_variants_from_file(input_path, file_name, variant_information, severe):
    """ reads in the variant from one patient file and adds the variants to the given dict. Key of dict is variant and value is 
        distribution over the categories.

    Args:
        input_path ( string ): path to input folder
        file_name ( string ): name of input file
        variant_information ( dict[string : array[integer]] ): dict with information about variants 
        severe ( boolean ): condition of patient
    """
    with open(input_path + "/" + file_name, 'r') as file:
        double_variants = set()
        for variant in file:
            
            # check for comments
            if variant[0] != '#':
                
                # get information of variant from line and create key
                content = variant.split('\t')
                chrom = content[0]
                pos = int(content[1])
                key = (chrom, pos)
                
                # avoid double counting of same variants
                double_variants.add(key)
                
                # add variant info to dict
                variant_information = add_variant_information(key, variant_information, severe, double_variants)
    return variant_information

def add_variant_information(key, variants, severe, double_variants):
    """ adds information about a single variant to the dict of variant information

    Args:
        key ( (string, integer) ): key of dict of a variant
        variants ( dict[(string, integer)] : array[integer] ): dict with variants information
        severe ( boolean ): condition of patient
        double_variants ( set[(string, integer)] ): already found variants 
    """
    # variant in dict (-> already found in other .vcf file)
    if key in variants and key not in double_variants:
        # update variant info
        value = variants.get(key)
        if severe:
            value[0] += 1
        else:
            value[2] += 1
    else:
        # save new variant
        if severe:
            variants[key] = [1, 0, 0, 0]
        else:
            variants[key] = [0, 0, 1, 0]
    return variants


print("--- START PROGRAM SUMMARYPROMOTERRESULTS.PY ---")

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

# print args
print(input_path)
print(info_file_path)
print(output_file_path)


# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS FOR PROMTER VCF FILES ---")

counter_pat, id_and_condition = get_patient_id(input_path, pattern)

# read conditions from Q001H_sample_preparations_20230803115337.tsv
# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = add_conditions(id_and_condition, info_file_path)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")
print("--- READ IN VARIANTS FROM PROMTER VCF FILES ---")

# dict entry: key = (chr, pos), value = (severePos, severeNeg, notSeverePos, not SevereNeg)
variant_information = {}
# iterate through all files from input folder
for file_name in os.listdir(input_path):
    print("Working in file: " + file_name)
    # get lab id and condition
    lab_id = extract_id(file_name, pattern)
    severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
    
    if  id_and_condition[lab_id] != "" and lab_id != "nullID":
        # add variants to dict
        variant_information = add_variants_from_file(input_path, file_name, variant_information, severe)
                        
print("--- SUCCESFUL ---")
print("--- SAVE FILE TO " + output_file_path + " ---")

# fill up dict
variant_information = fill_up_dict(variant_information, counter_severe, counter_not_severe)

# save to file
with open(output_file_path, 'w') as output_file:
    c = comments_file(output_file_path, input_path, counter_severe, counter_not_severe, len(variant_information))
    output_file.write(c)
    for variant in variant_information:
        s = variant_to_string(variant, variant_information[variant])
        output_file.write(s + "\n")
print("--- SUCCESSFUL -> PROGRAMM TERMINATES ---")