# Bachlor Bio Info
# author Martin Kahabka
# goal: Use the given names of the DE genes from "features.txt" to find the ENSG-codes
# the corresponding codes should be in the used gff file

# creates two file: 
# "EnsgOfDeGenes.txt" contains only the ENSG-Codes, "None" if no ENSG-Code is found
# "EnsgAndNameOfDeGenes.txt" contains the ENSG-Codes together with the gene_name 

import argparse

def extractInfo(s: str, feature: str) -> str:
    # extracts the information about a feature from a string: 
    # has to be format nameFeature ... "informationToExtract"

    index_feature = s.find(feature)
    if (index_feature == -1):
        # feature not found
        return ""
    else:
        # point in s after feature
        start_point = index_feature + len(feature)
        left_bound = s.find("\"", start_point) + 1
        right_bound = s.find("\"", left_bound)

        # return substring
        return s[left_bound:right_bound]

parser = argparse.ArgumentParser(description='Find ENSG-codes from gene names')

# optional arguemnts
parser.add_argument('-f', '--path_file', type=str, default='features.txt', help='Path to the file containing gene names')
parser.add_argument('-g', '--path_gff', type=str, default='../../Reference/Homo_sapiens.GRCh37.87.gtf', help='Path to the gff file')
parser.add_argument('-o', '--only_ensg', type=bool, default=True, help='Only ENSG-Codes are saved in file')

# Parse the arguments
args = parser.parse_args()

# Access the values of the arguments
file_path = args.path_file
gff_path = args.path_gff
only_ensg = args.only_ensg

print("-- read in DE genes --")

with open(file_path) as file:
    namesDeGenes = []
    for line in file:
        # extract names of differentially expressed genes
        data = line.split('\n')
        namesDeGenes.append(data[0])

print("-- read in DE genes -- FINISHED --")

# ENSG-Code of the differential expressed genes, starting element is "None"
print("-- read in gff_file and search for ENSG of genes --")

ensgForDe = ["None"] * len(namesDeGenes)
# read in lines of gff file
with open(gff_path) as gff:
    print("open file")
    for line in gff:
        data = line.split('\t')

        # check for commands and wrong formatted lines
        if line[0] != '#' and len(data) >= 9:
            # the additional info is saved an has to be extracted
            geneInfo = data[8]
            geneENSG = extractInfo(geneInfo, "gene_id")
            geneName = extractInfo(geneInfo, "gene_name")
            
            # check and possibly save ENSG
            if geneName in namesDeGenes:
                print("Found ENSG for gene: " + geneName)
                indexDeGene = namesDeGenes.index(geneName)
                ensgForDe[indexDeGene] = geneENSG

print("-- read in gff_file and search for ENSG of genes -- FINISHED --")

# print result
print("-- RESULT --")

for i in range(len(ensgForDe)):
    if ensgForDe[i] == "None":
        print("For gene " + namesDeGenes[i] + ": No ENSG found")
    else:
        print("For gene " + namesDeGenes[i] + ":  ENSG: " + ensgForDe[i])

# save info in files
print("-- SAVING IN FILE --")
filename_ensg = "EnsgOfDeGenes.txt"
filename_ensg_name  = "EnsgAndNameOfDeGenes.txt"


s_ensg = ""
s_ensg_name = ""

for i in range(len(ensgForDe)):
    # string for filename_ensg
    new_line_name = ensgForDe[i] + "\n"
    s_ensg = s_ensg + new_line_name

    # string for filename_ensg_name
    new_line_name = namesDeGenes[i] + "\t" + ensgForDe[i] + "\n"
    s_ensg_name = s_ensg_name + new_line_name

f = open(filename_ensg, "w")
f.write(s_ensg)
f.close

print("saved ensg codes under " + filename_ensg)

f = open(filename_ensg_name, "w")
f.write(s_ensg_name)
f.close

print("saved ensg with gene names under " + filename_ensg_name)