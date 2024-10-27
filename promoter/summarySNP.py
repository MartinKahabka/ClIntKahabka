# author: Martin Kahabka
import argparse
import os
import utils.summary_utils as s_utils

def fill_up_dict(dict_variant, num_severe, num_not_severe) -> dict:
    """ fills in the number of patients with no variants (posNegative and negNegativ) 
        based on the total number of patients and the number of patients with variants

    Args:
        dict_variant ( dict[string : array[integer]] ): dict with variant information
        num_severe (_type_): total number of patients with severe condition
        num_not_severe (_type_): total number of patients without severe condition

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
                
                # add variant info to dict
                variant_information = add_variant_information(key, variant_information, severe, double_variants)
                
                # avoid double counting of same variants
                double_variants.add(key)
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


print("--------------- START SUMMARYSNP.PY ---------------")

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
print("Input arguments")
print("Data files: " + input_path)
print("Patients info: " + info_file_path)
print("Output path: " + output_file_path)



# get IDs of patients
print("--- READ IN LAB IDS OF PATIENTS ---")

counter_pat, id_and_condition = s_utils.get_patient_id(input_path, pattern)

# read conditions from Q001H_sample_preparations_20230803115337.tsv
# extract lab_ID from patients
print("--- READ IN CONDITIONS OF PATIENTS ---")

id_and_condition, counter_severe, counter_not_severe = s_utils.add_conditions(id_and_condition, info_file_path, id_and_condition)

print("--- SUCCESFUL: number of severe/not severe patients: " + str(counter_severe) + "/" + str(counter_not_severe) + " ---")
print("--- READ IN VARIANTS FROM PROMTER VCF FILES ---")

# dict entry: key = (chr, pos), value = (severePos, severeNeg, notSeverePos, not SevereNeg)
variant_information = {}
# iterate through all files from input folder
for file_name in os.listdir(input_path):
    print("--- WORKING ON FILE: " + file_name + " ---")
    # get lab id and condition
    lab_id = s_utils.extract_id(file_name, pattern)
    severe = (True if id_and_condition[lab_id] == "COVID severe" else False)
    
    if  id_and_condition[lab_id] != "" and lab_id != "nullID":
        # add variants to dict
        variant_information = add_variants_from_file(input_path, file_name, variant_information, severe)
                        
# add number of patients where variants does not occur
variant_information = fill_up_dict(variant_information, counter_severe, counter_not_severe)

# save results to file
with open(output_file_path, 'w') as output_file:
    c = s_utils.comments_file(output_file_path, input_path, counter_severe, counter_not_severe, len(variant_information))
    output_file.write(c)
    for variant in variant_information:
        s = variant_to_string(variant, variant_information[variant])
        output_file.write(s + "\n")

print("--- SAVE FILE TO: " + output_file_path + " ---")
print("--------------- FINISHED SUMMARYSNP.PY ---------------")