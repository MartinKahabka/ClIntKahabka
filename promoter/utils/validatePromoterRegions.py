# author Martin Kahabka, 26.06.24
# input: 
# 1. file with promoter information in tsv format [chromsomeNumber, startPromoter, endPromoter, namePromoter, _, _]
# 2. gtf file of human genome (here GRCh37) in tsv format, for more infos see 'Homo_sapiens.GRCh37.87.gtf' file
#
# output:
# for every promoter determine if this promoter lies in the boundaries of the corresponding gene (according to the gtf file)
# either outputs:
# 1.    IN: + namePromoter
# 2.    OUT: + namePromoter
#       location relative to annotated gene region (up/downstream)
#       information of gene from gtf file
#       information to absolute location of the promoter (chr, start, end)
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


# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Script to validate promoter regions.')

# Add command line arguments
parser.add_argument('-g', '--gff_file', help='Path to the gtf file', default="Homo_sapiens.GRCh37.87.gtf")
parser.add_argument('-p', '--pro_pos_file', help='Path to the promoter file', default="GRCh37_promoterChrPos.bed")
parser.add_argument('-s', '--save_to_file', help='Should the output be saved to file', default=True)

# Parse the command line arguments
args = parser.parse_args()

# Access the values of the command line arguments
ggf_path = args.gff_file
promoter_reg_file = args.pro_pos_file
save_in_file = args.save_to_file

# read in file of promoter regions
promoter_info = {}
gene_to_promoter = {}
gene_names = set()
with open(promoter_reg_file, 'r') as file:
    for line in file:
        content = line.split('\t')
        chrNumber = int(content[0][3:])
        content[1] = int(content[1])
        content[2] = int(content[2])
        # add chromosome number and start/end site
        
        promoter_name = content[3]
        gene_name = content[3][:len(content[3])-2]
        
        # save promoter regions to gene name
        if gene_name in gene_to_promoter:
            gene_to_promoter[gene_name].append(promoter_name)
        else:
            gene_to_promoter[gene_name] = [promoter_name]
        
        # save promoter infos to promoter region name
        promoter_info[promoter_name] = (chrNumber, content[1], content[2])
        
        gene_names.add(gene_name)

# compare to gtf
report = True
str_to_file = ""
# open gtf file
with open(ggf_path, 'r') as file:
    for line in file:
        if line[0] not in ['#']:
            
            content = line.split('\t')
            try:
                content[0] = int(content[0])
            except:
                pass
            content[3] = int(content[3])
            content[4] = int(content[4])
            
            # get gene name of gft line
            gene_name = extractInfo(content[8], "gene_name")
            
            if gene_name in gene_names and content [2] == "gene":
                for promoter in gene_to_promoter.get(gene_name):
                    
                    info = promoter_info.get(promoter)
                    s  = ""
                    # check if same chromosome and if region of promotor lies in annotated gene region
                    if (info[0] == content[0]) and (info[1] >= content[3]) and (info[2] <= content[4]):
                        # everything is fine
                        s = "IN: " + str(promoter)
                        print(s)
                    else:
                        # promoter lies outside of annotated gene region
                        # print, saves string and poss. print to file
                        s = "OUT: " + str(promoter)
                        print(s)
                        # promoter lies downstream of annotated gene region
                        if info[1] < content[3]:
                            txt = "Downstream, diff: " + str(content[3] - info[1])
                            print(txt)
                            s += "\n" + txt
                        else:
                        # promoter lies downstream of annotated gene region
                            txt = "Upstream, diff: " + str(info[1] - content[4])
                            print(txt)
                            s += "\n" + txt
                        # print report
                        if report:
                            print(line)
                            print(info)
                            
                            s += "\n" + str(line) + "\n" + str(info)
                    
                    if save_in_file:
                        str_to_file += s + "\n"

with open("resultValidationPromoter.txt", 'w') as file:
    file.write(str_to_file)
                    