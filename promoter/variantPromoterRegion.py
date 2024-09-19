# author martin kahabka
# use: based on a list of variants and list of promoter. Every variant that is not in range of a promoter is filtered out.
#      The variants within a promoter are saved in a file
#      Also counts the number of variants per promoter and saves them in a file
import argparse
from functools import cmp_to_key
import os

class Region:
    def __init__(self, c, p):
        self.chrom = c
        self.pos = p
    
        
class Variant(Region):
    def __init__(self, c, p, l):
        super().__init__(c, p)
        self.line_content = l
    
    def identifier(self):
        """
        returns unique identifier for variant
        """
        return self.chrom + " " + str(self.pos)

class Promoter(Region):
    def __init__(self, c, p, l, n, num):
        super().__init__(c, p)
        self.line_content = l
        self.name = n
        self.num_variants = num
    
    def promToString(self):
        prom_as_array = [self.chrom, str(self.pos), self.name, str(self.num_variants)]
        return '\t'.join(prom_as_array)

"""
Converts a string identifier of a chromosome to and integer. 
Converts X, Y and MT string to numbers 23-25

@param chrom: string identifier of chromosome
@return: integer identifier of chromosome
"""

def chrToNumber(chrom : str):
    """ Converts a string identifier of a chromosome to and integer. \n
        Converts X, Y and MT string to numbers 23-25

    Args:
        chrom (string): string identifier of chromosome

    Returns:
        integer: integer identifier of chromosome
    """
    rest = chrom[3:]
    if rest == 'X':
        return 23
    if rest == 'Y':
        return 24
    if rest == "MT":
        return 25
    else:
        return int(rest)

def sortGenePos(x1 : Region, x2 : Region):
    """ Compares the position of two genome positions (chromosome and position) over the whole genome.

    Args:
        x1 (Region): first region to compare
        x2 (Region): second region to compar

    Returns:
        integer: 1 if x1 is bigger that x2 \n
                 -1 if x1 is lower that x2 \n
                 0 else
    """
    # get chr
    chr1 = chrToNumber(x1.chrom)
    chr2 = chrToNumber(x2.chrom)
    # lower chromosome, less big
    if chr1 < chr2:
        return -1
    elif chr1 > chr2:
        return 1
    # same chromosome
    else:
        pos1 = x1.pos
        pos2 = x2.pos
        if pos1 < pos2:
            return -1
        elif pos1 > pos2:
            return 1
    # pos are exactly the same
    return 0

def variantInBound(start : int, end : int, prom : Region, variant : Region) -> str:
    """ Compares the range of a given promoter to the position of a variant 
        and returns whether the variant is in, below or above the promoter

    Args:
        start (integer): lower bound of promoter (upstream)
        end (integer): upper bound of promoter (downstream)
        prom (Region): 
        variant (Region): _description_

    Returns:
        str: "bigger" if promoter range is higher that variant \n
             "smaller" if promoter range is lower that variant \n
             "in" if variant is in promoter range
    """
    # 1 if promoter is bigger that variant, -1 if otherwise, 0 else
    relPos = sortGenePos(prom, variant)
    promPos = prom.pos
    variantPos = variant.pos
    # check if chromosoms are equal
    if chrToNumber(prom.chrom) == chrToNumber(variant.chrom): 
        # define bound and check
        downBound = promPos - start
        upBound = promPos +  end
        if  downBound <= variantPos <= upBound:
            return "in"
        else:
            return "bigger" if relPos == 1 else "smaller"
    else:
        return "bigger" if relPos == 1 else "smaller"
    
def write_output_comment(path, n, vcf_p, promoter_p, start, end):
    """ Returns the comment for the output file

    Args:
        path (string): full output path to file
        n (string): name of process
        vcf_p (string): path of .vcf file of patient
        promoter_p (string): path to promoter file
        start (integer): range downstream of TSS
        end (integer): range upstream of TSS

    Returns:
        string: concatentation informations used as comment for output file
    """
    c = "# Name of file: " + path + "\n"
    c += "# Name of process: " + n + "\n"
    c += "# Patient vcf file path: " + vcf_p + "\n"
    c += "# Promoter region file path:" + promoter_p + "\n"
    c += "# Length Down/Upstream region of promoter TSS side: " + str(start) + "/" + str(end) + "\n"
    return c

def readInPromoter(prom_path: str):
    """ Read in file and returns list of promoters from this file

    Args:
        prom_path (string): Path to file with promoters 

    Returns:
        array[Promoter]: Array containing the promoters
    """
    with open(prom_path, 'r') as file:
        # get promoters regions of all promoters
        promoter_regions = []
        for line in file:
            if line[0] != '#':
                contain = line.split('\t')
                chrom = contain[0]
                pos = int(contain[1])
                name_prom = contain[3]
                num_of_vars = 0
                # add to promoter regions
                p = Promoter(chrom, pos, line, name_prom, num_of_vars)
                promoter_regions.append(p)
    return promoter_regions

print("--- START PROGRAMM VARIANTPROMOTERREGION.PY ---")

# arguments of running the code
parser = argparse.ArgumentParser(prog='variantPromoterRegion.py', description='Description of your script')
parser.add_argument('-n', '--name_process', help="one unique identifer of the process", required=False)
parser.add_argument('-o', '--output_dir',  help="where the output should be saved", required=False)
parser.add_argument('-v', '--vcf_path',  help="path to the vcf file of the patient", required=False)
parser.add_argument('-p', '--path_promoters', help="path to the file with the promoter regions", required=False)
parser.add_argument('-s', '--start', help="number of bases downstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-e', '--end', help="number of bases upstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-f', '--output_prom_path', help="output path for variant sum of promoter", required=False)

# parse arguments to vars
args = parser.parse_args()
name = args.name_process
output_path = args.output_dir
vcf_path = args.vcf_path
promoter_path = args.path_promoters
start_prom = int(args.start)
end_prom = int(args.end)
output_prom_path = args.output_prom_path

# read in promoter regions for file as promoter classes
print("--- READ IN PROMOTER FILE ---")
promoter_regions = readInPromoter(promoter_path)

            
# sort in case promoter regions aren't sorted
sorter = cmp_to_key(sortGenePos)
promoter_regions.sort(key = sorter)
print("--- SUCCESS ---")

# name output file
filename_vcf = os.path.basename(vcf_path)
full_output_path = os.path.join(output_path, name + "_promoterVcfs_" + filename_vcf)
full_sum_output_path = os.path.join(output_prom_path,  name + "_variantSum_" + filename_vcf)

print("--- START LOOKING FOR VCFS IN PROMOTER REGIONS IN: " + vcf_path + " ---")
# read in vcf of patient
with open(vcf_path, 'r') as vcf_file, open(full_output_path, 'w') as filter_vcfs_file:
    # write information into promoter vcf file
    filter_vcfs_file.write(write_output_comment(os.path.basename(full_output_path), name, vcf_path, promoter_path, start_prom, end_prom))
    
    # keep pointer on promoter, due to sorted arrays in O(n), n : num Vcfs in vcf file
    pointer_promoter = 0
    previous_variants = set()
    for line in vcf_file:
        if line[0] != '#' and line[0] != '\n':
            contain = line.split('\t')
            
            # define variant/promoter
            variant = Variant(contain[0], int(contain[1]), line)
            current_promoter = promoter_regions[pointer_promoter]
            
            # returns smaller, in or bigger. Relative pos of promoter to variant
            relPos = variantInBound(start_prom, end_prom, current_promoter, variant)
            
            # get current promoter
            while relPos == "smaller" and pointer_promoter+1 != len(promoter_regions):
                # variant is upstream/higher chrom that promoter, step to next
                pointer_promoter += 1
                relPos = variantInBound(start_prom, end_prom, promoter_regions[pointer_promoter], variant)
                
            # update promoter
            current_promoter = promoter_regions[pointer_promoter]
        
            if relPos == "in" and variant.identifier() not in previous_variants:
                print("found vcf in promoter " + current_promoter.name + " on " + variant.chrom + ", pos: " + str(variant.pos))
                # write to file
                filter_vcfs_file.write(variant.line_content)
                # update counter and add to known variants to void copies
                previous_variants.add(variant.identifier())
                current_promoter.num_variants += 1
                
# save sum of variants per promoter to file
with open(full_sum_output_path, "w") as promoter_sum_file:
    for promoter in promoter_regions:
        promoter_sum_file.write(promoter.promToString() + "\n")

print("SUCCESS: PROGRAMM FINISHED")