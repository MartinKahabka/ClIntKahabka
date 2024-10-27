import re
import os

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

def add_conditions(dict_id, filename, id_and_condition):
    """ matches conditions (severe/not severe) to the the patients ids

    Args:
        dict_id ( dict[string : string] ): input dict without conditions
        filename (string): path to file that contains the condition of the patients
        id_and_condition: dict[string : string] dict without conditions 

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