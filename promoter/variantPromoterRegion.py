import argparse

parser = argparse.ArgumentParser(prog='variantPromoterRegion.py', description='Description of your script')

parser.add_argument('-n', '--name_process', default='name_process', help="one unique identifer of the process")
parser.add_argument('-o', '--output_dir', default='output_dir', help="where the output should be saved")
parser.add_argument('-v', '--vcf_path', default='vcf_path', help="path to the vcf file of the patient")
parser.add_argument('-p', '--path_promoters', default='path_promoters', help="path to the file with the promoter regions")
parser.add_argument('-s', '--start', default='start', help="number of bases downstream of the promoter TSS that are considered promoter region")
parser.add_argument('-e', '--end', default='end', help="number of bases upstream of the promoter TSS that are considered promoter region")

args = parser.parse_args()

name = args.n
output_path = args.o
vcf_path = args.v
promoter_path = args.p
start_prom = args.s
end_prom = args.e


# safe promoter regions as tuple (chromosome, position TSS)
promoter_regions = []

with open(promoter_path, 'r') as file:
    # get promoters regions of all promoters
    for line in file:
        if line[0] != '#':
            contain = line.split('\t')
            chrom = contain[0]
            pos = contain[1]
            # add to promoter regions
            promoter_regions.add(chrom, pos)