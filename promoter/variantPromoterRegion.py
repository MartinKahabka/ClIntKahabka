import argparse
from functools import cmp_to_key

# run with comment
#python3

def chrToNumber(chrom):
    if chrom == 'X':
        return 23
    if chrom == 'Y':
        return 24
    else:
        return int(chrom[3:])

def sortPromoter(x1, x2):
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
start_prom = args.start
end_prom = args.end


# safe promoter regions as tuple (chromosome, position TSS)
promoter_regions = []

with open(promoter_path, 'r') as file:
    # get promoters regions of all promoters
    for line in file:
        if line[0] != '#':
            contain = line.split('\t')
            chrom = contain[0]
            pos = int(contain[1])
            # add to promoter regions
            promoter_regions.append((chrom, pos))
            
# if regions arent sorted
sorter = cmp_to_key(sortPromoter)
promoter_regions.sort(key = sorter)
