import argparse
from functools import cmp_to_key
import os

# run with comment
#python3
class vcf_file:
    def __init__(self, input_path):
        self.file = open(input_path, "r")
        self.size = os.path.getsize(input_path)
        self.end_byte = self.calculate_end_byte()
        self.start_byte = self.calculate_start_byte()
        
    def calculate_start_byte(self):
        self.file.seek(0)
        while self.file.read(1) == '#':
            self.file.readline()
        return self.file.tell()
    
    def calculate_end_byte(self):
        offset = 1
        self.file.seek(self.size - offset)
        while self.file.read(1) == '\n':
            offset += 1
            self.file.seek(self.size - offset)
        return self.size - offset
    
    def getLine(self, byte):
        self.file.seek(byte, 0)
        # at linebreak return file
        if self.file.read(1) != "\n":
            # go back to last linebreak
            self.file.seek(byte, 0)
            start_offset = 0
            while self.file.read(1) != "\n":
                start_offset += 1
                self.file.seek(byte - start_offset)
        start_index = self.file.tell()
        # create vars for line class
        line = self.file.readline()
        start_index = start_index
        end_index = start_index + len(line)
        return vcf_line(line, start_index, end_index)
        
    def nextLine(self, line_of_vcf):
        byte = line_of_vcf.end_byte + 1
        return self.getLine(byte)
        
    def previousLine(self, line_of_vcf):
        byte = line_of_vcf.start_byte - 2
        return self.getLine(byte)
    
class output_file:
    def __init__(self, path):
        self.path = path
        self.file = open(path, "w")
        
    def write(self, s):
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
        # get chr
        chr1 = chrToNumber(self.chr)
        chr2 = chrToNumber(second_gene[0])
        # lower chromosome, less big
        if chr1 < chr2:
            return -1
        elif chr1 > chr2:
            return 1
        # same chromosome
        else:
            pos1 = self.pos
            pos2 = second_gene[1]
            if pos1 < pos2:
                return -1
            elif pos1 > pos2:
                return 1
        # pos are exactly the same
        return 0
    
    def lineInPromoter(self, start, end, prom):
        # -1 if promoter is downstream of variant, 1 if otherwise, 0 else
        relPos = self.compareGenePos(prom)
        promPos = prom[1]
        variantPos = self.gene[1]
        # check if chromosoms are equal
        if chrToNumber(prom[0]) == chrToNumber(self.gene[0]): 
            # define bound and check
            downBound = promPos - start
            upBound = promPos +  end
            if  downBound <= variantPos <= upBound:
                return "in"
        return "out"

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
                
def getCommand(output_path, name, vcf_path, promoter_path, start, end) -> str:
    command = ""
    command += "# Name of file: " + os.path.basename(output_path) + "\n"
    command += "# Name of process: " + name + "\n"
    command += "# Patient vcf file path: " + vcf_path + "\n"
    command += "# Promoter region file path:" + promoter_path + "\n"
    command += "# Length Down/Upstream region of promoter TSS side: " + str(start) + "/" + str(end) + "\n"
    return command

def fastAnalyis_refac(input_file, result_file, prom, upstream, downstream):
    # binary search
    start = input_file.start_byte
    end = input_file.end_byte
    while start < end:
        index = start + int((end - start) / 2)
        line = input_file.getLine(index)

        if line.isGene:
            rel_pos = line.compareGenePos(prom)
            
            #print(str(start) + "||" + str(end))
            #print(line.raw_data)
            #print(prom)
            #print(rel_pos)
            # gene has bigger position
            if rel_pos == -1:
                start = index + 1
            # promoter has bigger position
            elif rel_pos == 1:
                end = index - 1
            # position of gene and promoter is identical
            else:
                index = start
                break

    # search lowest line in promoter
    lowest_line = line
    previous_line = input_file.previousLine(lowest_line)
    # check if line lies in boundries of promoter and file
    while previous_line.lineInPromoter(upstream, downstream, prom) == "in" and previous_line.start_byte >= input_file.start_byte:        
        previous_line = input_file.previousLine(lowest_line)
        lowest_line = previous_line
    lowest_line = input_file.nextLine(lowest_line)
    
    # find all lines that lie in promoter
    while lowest_line.isGene and lowest_line.lineInPromoter(upstream, downstream, prom) == "in":
        
        result_file.write(lowest_line.raw_data)
        lowest_line = input_file.nextLine(lowest_line)
        
    
    print("done")

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

print("--- START PROGRAMM VARIANTPROMOTERREGION_REFACT.PY ---")

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
            # add to promoter regions
            promoter_regions.append((chrom, pos, name_prom))
            
# sort in case promoter regions aren't sorted
sorter = cmp_to_key(sortGenePos)
promoter_regions.sort(key = sorter)
print("--- SUCCESS ---")

# name output file
filename_vcf = os.path.basename(vcf_path)
full_output_path = os.path.join(output_path, name + "_promoterVcfs_" + filename_vcf)

print("--- START LOOKING FOR VCFS IN PROMOTER REGIONS IN: " + vcf_path + " ---")
# read in vcf of patient
input_file = vcf_file(vcf_path)
result_file = output_file(full_output_path)
# write command of output file
result_file.write("# Name of file: " + os.path.basename(full_output_path) + "\n")
result_file.write("# Name of process: " + name + "\n")
result_file.write("# Patient vcf file path: " + vcf_path + "\n")
result_file.write("# Promoter region file path:" + promoter_path + "\n")
result_file.write("# Length Down/Upstream region of promoter TSS side: " + str(start_prom) + "/" + str(end_prom) + "\n")

for promoter in promoter_regions:
    fastAnalyis_refac(input_file, result_file, promoter, start_prom, end_prom)


print("SUCCESS: PROGRAMM FINISHED")