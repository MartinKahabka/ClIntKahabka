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

class ROI_region(Region):
    def __init__(self, c, p, l, n, num):
        super().__init__(c, p)
        self.line_content = l
        self.name = n
        self.num_variants = num
    
    def promToString(self):
        prom_as_array = [self.chrom, str(self.pos), self.name, str(self.num_variants)]
        return '\t'.join(prom_as_array)

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

def readInROIs(prom_path: str):
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
                p = ROI_region(chrom, pos, line, name_prom, num_of_vars)
                promoter_regions.append(p)
    return promoter_regions

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