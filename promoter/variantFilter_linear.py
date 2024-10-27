# author martin kahabka
# use: based on a list of variants and list of promoter. Every variant that is not in range of a promoter is filtered out.
#      The variants within a promoter are saved in a file. Uses linear algorithm.
#      Also counts the number of variants per promoter and saves them in a file
import argparse
from functools import cmp_to_key
import os
from utils.variantFilter_utils import Region, Variant
import utils.variantFilter_utils as vF_utils

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
    relPos = vF_utils.sortGenePos(prom, variant)
    promPos = prom.pos
    variantPos = variant.pos
    # check if chromosoms are equal
    if vF_utils.chrToNumber(prom.chrom) == vF_utils.chrToNumber(variant.chrom): 
        # define bound and check
        downBound = promPos - start
        upBound = promPos +  end
        if  downBound <= variantPos <= upBound:
            return "in"
        else:
            return "bigger" if relPos == 1 else "smaller"
    else:
        return "bigger" if relPos == 1 else "smaller"

print("--------------- START VARIANTFILTER_LINEAR.PY ---------------")
# read and parse input parameters
parser = argparse.ArgumentParser(prog='variantPromoterRegion.py', description='Description of your script')
parser.add_argument('-n', '--name_process', help="one unique identifer of the process", required=False)
parser.add_argument('-o', '--output_dir',  help="where the output should be saved", required=False)
parser.add_argument('-v', '--vcf_path',  help="path to the vcf file of the patient", required=False)
parser.add_argument('-p', '--path_promoters', help="path to the file with the promoter regions", required=False)
parser.add_argument('-s', '--start', help="number of bases downstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-e', '--end', help="number of bases upstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-f', '--output_prom_path', help="output path for sum of variants per region", required=False)

args = parser.parse_args()
name = args.name_process
output_path = args.output_dir
vcf_path = args.vcf_path
ROI_path = args.path_promoters
start_region = int(args.start)
end_region = int(args.end)
output_sum_path = args.output_prom_path

# print args
print("Input arguments")
print("name: " + name)
print("vcf file: " + vcf_path)
print("regions of interest: " + ROI_path)
print("boundary upstream: " + str(start_region))
print("boundary downstream: " + str(end_region))
print("Output path filtered vcfs: " + output_path)
print("Output path sum of variants per region: " + output_sum_path)

# read in regions of interest from file as ROI_region classes
print("--- READ IN PROMOTER FILE ---")
regions_of_interest = vF_utils.readInROIs(ROI_path)

            
# sort in case ROIs regions aren't sorted
sorter = cmp_to_key(vF_utils.sortGenePos)
regions_of_interest.sort(key = sorter)

# create names of output files
filename_vcf = os.path.basename(vcf_path)
full_output_path = os.path.join(output_path, name + "_promoterVcfs_" + filename_vcf)
full_sum_output_path = os.path.join(output_sum_path,  name + "_variantSum_" + filename_vcf)

print("--- START LOOKING FOR VCFS IN PROMOTER REGIONS IN: " + vcf_path + " ---")
pointer_regions = 0
previous_variants = set()

# read in vcf of patient
with open(vcf_path, 'r') as vcf_file, open(full_output_path, 'w') as filter_vcfs_file:
    # write information into promoter vcf file
    filter_vcfs_file.write(vF_utils.write_output_comment(os.path.basename(full_output_path), name, vcf_path, ROI_path, start_region, end_region))
    
    # iterate over variants
    for line in vcf_file:
        if line[0] != '#' and line[0] != '\n':
            contain = line.split('\t')

            variant = Variant(contain[0], int(contain[1]), line)
            current_region = regions_of_interest[pointer_regions]
            
            # returns smaller, in or bigger. Relative pos of promoter to variant
            relPos = variantInBound(start_region, end_region, current_region, variant)
            
            # get current promoter
            while relPos == "smaller" and pointer_regions+1 != len(regions_of_interest):
                # variant is upstream/higher chrom that promoter, step to next
                pointer_regions += 1
                relPos = variantInBound(start_region, end_region, regions_of_interest[pointer_regions], variant)
                
            # update promoter
            current_region = regions_of_interest[pointer_regions]
        
            if relPos == "in" and variant.identifier() not in previous_variants:
                # write to file
                filter_vcfs_file.write(variant.line_content)
                # update counter and add to known variants to void copies
                previous_variants.add(variant.identifier())
                current_region.num_variants += 1
                
# save sum of variants per promoter to file
with open(full_sum_output_path, "w") as promoter_sum_file:
    for promoter in regions_of_interest:
        promoter_sum_file.write(promoter.promToString() + "\n")

print("--------------- FINISHED VARIANTFILTER_LINEAR.PY ---------------")