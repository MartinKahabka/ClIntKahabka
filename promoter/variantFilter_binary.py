# author martin kahabka
# use: based on a list of variants and list of promoter. Every variant that is not in range of a promoter is filtered out.
#      The variants within a promoter are saved in a file. Uses binary search algorithm.
#      Also counts the number of variants per promoter and saves them in a file
import argparse
from functools import cmp_to_key
import os

class vcf_file:
    def __init__(self, input_path):
        self.file = open(input_path, "r")
        self.size = os.path.getsize(input_path)
        self.end_byte = self.calculate_end_byte()
        self.start_byte = self.calculate_start_byte()
        
    def calculate_start_byte(self):
        """ calculate the first byte of the file excluding bytes that are part of the comment 

        Returns:
            integer: first byte that is not part of the comment
        """
        self.file.seek(0)
        while self.file.read(1) == '#':
            self.file.readline()
        
        return self.file.tell() - 2
    
    def calculate_end_byte(self):
        """ calculate the last byte of the file excluding byte that  are part of overhanging newline chars

        Returns:
            integer: last byte that is not a newline char
        """
        offset = 1
        self.file.seek(self.size - offset)
        while self.file.read(1) == '\n':
            offset += 1
            self.file.seek(self.size - offset)
        return self.size - offset
    
    def getLine(self, byte):
        """ based on a byte from the file returns a line in the file, defined by two (start and end) newline characters

        Args:
            byte (integer): byte that is part of the line

        Returns:
            string: found line in the file
        """
        self.file.seek(byte, 0)
        # at linebreak return file
        if self.file.read(1) != "\n":
            # go back to last linebreak
            self.file.seek(byte, 0)
            start_offset = 0
            while self.file.read(1) != "\n" and start_offset < byte:
                start_offset += 1
                self.file.seek(byte - start_offset)
        start_index = self.file.tell()
        # create vars for line class
        line = self.file.readline()
        start_index = start_index
        end_index = start_index + len(line)
        return vcf_line(line, start_index, end_index)
        
    def nextLine(self, line_of_vcf):
        """ given a line in the file, returns the next line from the file

        Args:
            line_of_vcf (string): line in file

        Returns:
            string: next line in file
        """
        byte = line_of_vcf.end_byte + 1
        return self.getLine(byte)
        
    def previousLine(self, line_of_vcf):
        """ given a line in the file, returns the previous line from the file

        Args:
            line_of_vcf (string): line in file

        Returns:
            string: previous line of file
        """
        byte = line_of_vcf.start_byte - 2
        return self.getLine(byte)
        
class output_file:
    def __init__(self, path):
        self.path = path
        self.file = open(path, "w")
        
    def write(self, s):
        """ writes to output file

        Args:
            s (string): string to write to file
        """
        self.file.write(s)

class vcf_line:
    def __init__(self, contains, start_byte, end_byte):
        self.raw_data = contains
        self.array_data = contains.split("\t")
        self.start_byte = start_byte
        self.end_byte = end_byte
        try:
            self.chr = self.array_data[0]
            self.pos = int(self.array_data[1])
            self.gene = (self.array_data[0], int(self.array_data[1]))
            self.isGene = True
        except:
            self.gene = None
            self.chr = None
            self.pos = None
            self.isGene = False
            
    def getStringLine(self):
        return self.raw_data
    
    def getGenePosition(self):
        return self.array_data[0], int(self.array_data[1])
    
    def compareGenePos(self, second_gene):
        """ compares the position of the line to a given gene region

        Args:
            second_gene (Region): gene region from position comparison

        Returns:
            integer: 1 if line is bigger that region \n
                 -1 if line is lower that region \n
                 0 else
        """
        # get chr
        chr1 = chrToNumber(self.chr)
        chr2 = chrToNumber(second_gene.chrom)
        # lower chromosome, less big
        if chr1 < chr2:
            return -1
        elif chr1 > chr2:
            return 1
        # same chromosome
        else:
            pos1 = self.pos
            pos2 = second_gene.pos
            if pos1 < pos2:
                return -1
            elif pos1 > pos2:
                return 1
        # pos are exactly the same
        return 0
    
    def lineInPromoter(self, start, end, prom):
        """ checks whether this line is in the range of a promoter

        Args:
            start (integer): lower bound of promoter
            end (integer): upper bound of promoter
            prom (region): promoter to check for

        Returns:
            string: "in" if line is in promoter
                    "out" else
        """
        # -1 if promoter is downstream of variant, 1 if otherwise, 0 else
        promPos = prom.pos
        variantPos = self.pos
        # check if chromosoms are equal
        if chrToNumber(prom.chrom) == chrToNumber(self.chr): 
            # define bound and check
            downBound = promPos - start
            upBound = promPos +  end
            if  downBound <= variantPos <= upBound:
                return "in"
        return "out"

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

def chrToNumber(chrom):
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
                    
def getCommand(output_path, name, vcf_path, promoter_path, start, end) -> str:
    """ Returns the comment for the output file

    Args:
        output_path (string): full output path to file
        name (string): name of process
        vcf_path (string): path of .vcf file of patient
        promoter_path (string): path to promoter file
        start (integer): range downstream of TSS
        end (integer): range upstream of TSS

    Returns:
        string: concatentation informations used as comment for output file
    """
    command = ""
    command += "# Name of file: " + os.path.basename(output_path) + "\n"
    command += "# Name of process: " + name + "\n"
    command += "# Patient vcf file path: " + vcf_path + "\n"
    command += "# Promoter region file path:" + promoter_path + "\n"
    command += "# Length Down/Upstream region of promoter TSS side: " + str(start) + "/" + str(end) + "\n"
    return command

def binarySearch_line(start_byte, end_byte, variant_file, current_promoter):
    """ conducts a binary search to find the variant that is, in relation to the other variants, 
        the closes (above or below) the promoter

    Args:
        start_byte (integer): start byte of file content (excluding comment at beginning)
        end_byte (integer): end byte of file content (excluding newlines at end)
        variant_file (vcf_file): vcf file with variants
        current_promoter (promoter): promoter to which position to search for

    Returns:
        vcf_line: line closest (below or above) to the promoter
    """
    # start binary search
    while start_byte < end_byte:
        # get current variant by index
        index = start_byte + int((end_byte - start_byte) / 2)
        current_line = variant_file.getLine(index)

        # adjust search range
        if current_line.isGene:
            rel_pos = current_line.compareGenePos(current_promoter)
            
            # gene has bigger position
            if rel_pos == -1:
                start_byte = index + 1
            # promoter has bigger position
            elif rel_pos == 1:
                end_byte = index - 1
            # position of gene and promoter is identical
            else:
                index = start_byte
                return current_line
    return current_line

def find_lowest_variant(lowest_variant, variant_file, current_promoter, upper_bound, lower_bound):
    """ given a line finds the lowest variant in the file that still is in bound of the promoter

    Args:
        lowest_variant (vcf_line): variant to search below
        variant_file (vcf_file): vcf file with variants
        current_promoter (region): promoter with bounds
        upper_bound (integer): lower bound of promoter
        lower_bound (integer): upper bound of promoter

    Returns:
        vcf_line: lowest variant in bound of promoter or input variant if no variant below is in range
    """
    previous_variant = variant_file.previousLine(lowest_variant)
    
    # get most upstream variant in boundary of promoter
    while previous_variant.start_byte >= variant_file.start_byte and previous_variant.lineInPromoter(upper_bound, lower_bound, current_promoter) == "in":        
        lowest_variant = previous_variant
        previous_variant = variant_file.previousLine(lowest_variant)
        
    return lowest_variant

def find_variants_in_bound(lowest_variant, variant_file, current_promoter, upper_bound, lower_bound, previous_variants):
    """ given a variant finds all variants above this variant that are in bound of the promoter

    Args:
        lowest_variant (vcf_line): variant to search above
        variant_file (vcf_file): vcf file with variants
        current_promoter (promoter): promoter with bounds
        upper_bound (integer): lower bound of promoter
        lower_bound (integer): upper bound of promoter
        previous_variants (set[promoter]): previous variants that are already considered

    Returns:
        (set[promoter], array[vcf_line]): found variants that are in bound and updated previous variants
    """
    variants_in_bound = []
    
    # check for edge case
    if lowest_variant.lineInPromoter(upper_bound, lower_bound, current_promoter) != "in":
        lowest_variant = variant_file.nextLine(lowest_variant)
    
    while lowest_variant.isGene and lowest_variant.lineInPromoter(upper_bound, lower_bound, current_promoter) == "in":
        if lowest_variant.getGenePosition() not in previous_variants:
            variants_in_bound.append(lowest_variant)
            previous_variants.add(lowest_variant.getGenePosition())
            current_promoter.num_variants += 1
        lowest_variant = variant_file.nextLine(lowest_variant)
        
    return (previous_variants, variants_in_bound)

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

def variantFiltering_binary(input_file, result_file, prom_sum_file, prom, upstream, downstream, prev_vars):
    """ conducts a binary search for a given promoter. Save the found variants to a output file
        and the number of variant per promoter to another output file

    Args:
        input_file (vcf_file): vcf file with variants
        result_file (string): output file path for filtered variants
        prom_sum_file (string): output file path for sum of variants per promoter
        prom (promoter): promoter to search for variants
        upstream (integer): lower bound of promoter
        downstream (integer): upper bound of promoter
        prev_vars (set[promoter]): previous considered variants

    Returns:
        set[promoter]: updated previous considered variants
    """
    # binary search
    start = input_file.start_byte
    end = input_file.end_byte
    
    # check if vcf file is not empty
    if start < end:
        
        # closest line at promoter (below or above)
        lowest_line = binarySearch_line(start, end, input_file, prom)
        
        # get most upstream variant in boundary of promoter
        lowest_line = find_lowest_variant(lowest_line, input_file, prom, upstream, downstream)
        
        # search upwards for other variants and count sum of variants
        prev_vars, variants_in_bound = find_variants_in_bound(lowest_line, input_file, prom, upstream, downstream, prev_vars)
        
        # write found variants to output file
        for variant in variants_in_bound:
            result_file.write(variant.raw_data)
        
    # write amount of variant in promoter to output file
    prom_sum_file.write(prom.promToString() + '\n')
    
    return prev_vars

print("--- START PROGRAMM VARIANTPROMOTERREGION_REFACT.PY ---")

parser = argparse.ArgumentParser(prog='variantPromoterRegion.py', description='Description of your script')

parser.add_argument('-n', '--name_process', help="one unique identifer of the process", required=False)
parser.add_argument('-o', '--output_dir',  help="where the output should be saved", required=False)
parser.add_argument('-v', '--vcf_path',  help="path to the vcf file of the patient", required=False)
parser.add_argument('-p', '--path_promoters', help="path to the file with the promoter regions", required=False)
parser.add_argument('-s', '--start', help="number of bases downstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-e', '--end', help="number of bases upstream of the promoter TSS that are considered promoter region", required=False)
parser.add_argument('-f', '--output_prom_path', help="output path for variant sum of promoter", required=False)

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
full_vcf_output_path = os.path.join(output_path, name + "_promoterVcfs_" + filename_vcf)
full_sum_output_path = os.path.join(output_prom_path,  name + "_variantSum_" + filename_vcf)

print("--- START LOOKING FOR VCFS IN PROMOTER REGIONS IN: " + vcf_path + " ---")
# read in vcf of patient
input_file = vcf_file(vcf_path)
filtered_vcf_file = output_file(full_vcf_output_path)
promoter_sum_file = output_file(full_sum_output_path)

# write command of output file
filtered_vcf_file.write(write_output_comment(os.path.basename(full_vcf_output_path), name, vcf_path, promoter_path, start_prom, end_prom))

# to avoid doubles
previous_variants = set()
for promoter in promoter_regions:
    variantFiltering_binary(input_file, filtered_vcf_file, promoter_sum_file, promoter, start_prom, end_prom, previous_variants)

print("FILE SAVED TO: " + filtered_vcf_file.path)
print("SUCCESS: PROGRAMM FINISHED")