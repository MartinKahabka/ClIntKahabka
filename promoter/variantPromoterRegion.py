import argparse
from functools import cmp_to_key
import os

# run with comment
#python3

def chrToNumber(chrom):
    rest = chrom[3:]
    if rest == 'X':
        return 23
    if rest == 'Y':
        return 24
    if rest == "MT":
        return 25
    else:
        return int(rest)

def sortGenePos(x1, x2):
    # get chr
    chr1 = chrToNumber(x1[0])
    chr2 = chrToNumber(x2[0])
    # lower chromosome, less big
    if chr1 < chr2:
        return -1
    elif chr1 > chr2:
        return 1
    # same chromosome
    else:
        pos1 = x1[1]
        pos2 = x2[1]
        if pos1 < pos2:
            return -1
        elif pos1 > pos2:
            return 1
    # pos are exactly the same
    return 0

def variantInBound(start, end, prom, variant) -> str:
    # -1 if promoter is downstream of variant, 1 if otherwise, 0 else
    relPos = sortGenePos(prom, variant)
    promPos = prom[1]
    variantPos = variant[1]
    # check if chromosoms are equal
    if chrToNumber(prom[0]) == chrToNumber(variant[0]): 
        # define bound and check
        downBound = promPos - start
        upBound = promPos +  end
        if  downBound <= variantPos <= upBound:
            return "in"
        else:
            return "bigger" if relPos == 1 else "smaller"
    else:
        return "bigger" if relPos == 1 else "smaller"
    
# prom of type [chrom, int(pos), name_prom, int(num_of_vars)]
def promToString(prom):
    prom_as_string = [str(i) for i in prom]
    return '\t'.join(prom_as_string)

print("--- START PROGRAMM VARIANTPROMOTERREGION.PY ---")

parser = argparse.ArgumentParser(prog='variantPromoterRegion.py', description='Description of your script')

parser.add_argument('-n', '--name_process', help="one unique identifer of the process", required=False)
parser.add_argument('-o', '--output_dir',  help="where the output should be saved", required=False)
parser.add_argument('-v', '--vcf_path',  help="path to the vcf file of the patient", required=False)
parser.add_argument('-p', '--path_promoters', help="path to the file with the promoter regions", required=False)
parser.add_argument('-s', '--start', help="number of bases downstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-e', '--end', help="number of bases upstream of the promoter TSS that are considered promoter region", required=False)

args = parser.parse_args()

name = args.name_process
output_path = args.output_dir
vcf_path = args.vcf_path
promoter_path = args.path_promoters
start_prom = int(args.start)
end_prom = int(args.end)


# safe promoter regions as tuple (chromosome, position TSS)
promoter_regions = []

# read in promoter regions
print("--- READ IN PROMOTER FILE ---")
with open(promoter_path, 'r') as file:
    # get promoters regions of all promoters
    for line in file:
        if line[0] != '#':
            contain = line.split('\t')
            chrom = contain[0]
            pos = int(contain[1])
            name_prom = contain[3]
            num_of_vars = 0
            # add to promoter regions
            promoter_regions.append([chrom, pos, name_prom, num_of_vars])
            
# sort in case promoter regions aren't sorted
sorter = cmp_to_key(sortGenePos)
promoter_regions.sort(key = sorter)
print("--- SUCCESS ---")

# name output file
filename_vcf = os.path.basename(vcf_path)
full_output_path = os.path.join(output_path, name + "_promoterVcfs_" + filename_vcf)

print("--- START LOOKING FOR VCFS IN PROMOTER REGIONS IN: " + vcf_path + " ---")
# read in vcf of patient
with open(vcf_path, 'r') as vcf_file, open(full_output_path, 'w') as filter_vcfs_file:
    # write information into promoter vcf file
    filter_vcfs_file.write("# Name of file: " + os.path.basename(full_output_path) + "\n")
    filter_vcfs_file.write("# Name of process: " + name + "\n")
    filter_vcfs_file.write("# Patient vcf file path: " + vcf_path + "\n")
    filter_vcfs_file.write("# Promoter region file path:" + promoter_path + "\n")
    filter_vcfs_file.write("# Length Down/Upstream region of promoter TSS side: " + str(start_prom) + "/" + str(end_prom) + "\n")
    # keep pointer on promoter, due to sorted arrays in O(n), n : num Vcfs in vcf file
    pointer_promoter = 0
    previous_variants = set()
    for line in vcf_file:
        if line[0] != '#' and line[0] != '\n':
            contain = line.split('\t')
            # (chromosome, position, ref, alt)
            variant = (contain[0], int(contain[1]), contain[3], contain[4])
            # returns smaller, in or bigger. Relative pos of promoter to variant
            relPos = variantInBound(start_prom, end_prom, promoter_regions[pointer_promoter], variant)
            
            # step to next promoters if current is downstream
            while relPos == "smaller" and pointer_promoter+1 != len(promoter_regions):
                pointer_promoter += 1
                relPos = variantInBound(start_prom, end_prom, promoter_regions[pointer_promoter], variant)
        
            if relPos == "in" and variant not in previous_variants:
                print("found vcf in promoter " + promoter_regions[pointer_promoter][2] + " on " + promoter_regions[pointer_promoter][0] + ", pos: " + str(contain[1]))
                previous_variants.add(variant)
                promoter_regions[pointer_promoter][3] += 1
                filter_vcfs_file.write(line)
                
# save sum of variants per promoter to file
output_prom_path = "informationAndData/output_promoter/validate/variants_per_promoter"
full_sum_output_path = os.path.join(output_prom_path,  name + "_variantSum_" + filename_vcf)

with open(full_sum_output_path, "w") as promoter_sum_file:
    for promoter in promoter_regions:
        promoter_sum_file.write(promToString(promoter) +  '\n')

print("SUCCESS: PROGRAMM FINISHED")