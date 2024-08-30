# author martin kahabka
# use filters regions of interest from database (.bed) based on selection of gene names
import argparse

# input: name of region (promoter/enhancer/silencer name)
# output: name of gene corresponding to region (name_number -> name)
def regionToGeneName(region: str):
    name_gene = region.split("_")
    return name_gene[0]

print("--------------- START REGIONOFINTERESTFILTER.PY ----------------")

# input arguments of file
parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-d", "--database_path", help="Path to the database")
parser.add_argument("-g", "--gene_names", help="File with gene names")
parser.add_argument("-o", "--output_dir", help="Path of output file")
args = parser.parse_args()

# read in parameters
path_dataset_file = args.database_path
path_genes_file = args.gene_names
path_output_file = args.output_dir

# read gene names from file
print("read in gene names")

gene_names = set()
with open(path_genes_file, "r") as gene_names_file:
    for name in gene_names_file:
        name = name.strip()
        gene_names.add(name)

# read in dataset and filter
print("read in dataset and write to output file")
num_filtered = 0
with open(path_dataset_file, "r") as dataset, open(path_output_file, "w") as output_file:
    for region in dataset:
        # read in content of line
        region = region.strip()
        content = region.split("\t")
        
        # get name of gene from name of region ("gene_num" -> "gene")
        gene_name = regionToGeneName(content[3])
        
        # filter for regions of interest
        if gene_name in gene_names:
            output_file.write(region + "\n")
            num_filtered += 1

print("Output saved under " + path_output_file)
print("Amount of regions of interest: " +  str(num_filtered))
print("--------------- FINISHED REGIONOFINTERESTFILTER.PY ----------------")