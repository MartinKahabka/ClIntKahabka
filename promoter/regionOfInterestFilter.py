# author martin kahabka
import argparse

def regionToGeneName(region: str):
    name_gene = region.split("_")
    return name_gene[0]

# REPLACE with input parameter options
print("--------------- START REGIONOFINTERESTFILTER.PY ----------------")

parser = argparse.ArgumentParser(description="Process input directory")
parser.add_argument("-d", "--database_path", help="Path to the database")
parser.add_argument("-g", "--gene_names", help="File with gene names")
parser.add_argument("-o", "--output_dir", help="Path of output file")
args = parser.parse_args()

# read in parameters
path_dataset_file = args.database_path
path_genes_file = args.gene_names
path_output_file = args.output_dir

# read gene names into set
gene_names = set()

print("read in gene names")
with open(path_genes_file, "r") as gene_names_file:
    for name in gene_names_file:
        name = name.strip()
        gene_names.add(name)

print("read in dataset and write to output file")
counter_output = 0
with open(path_dataset_file, "r") as dataset, open(path_output_file, "w") as output_file:
    for line in dataset:
        line = line.strip()
        content = line.split("\t")
        
        # get name of gene from name of region ("gene_num" -> "gene")
        gene_name = regionToGeneName(content[3])
        
        # filter for regions of interest
        if gene_name in gene_names:
            output_file.write(line + "\n")
            counter_output += 1
            
print("Output saved under " + path_output_file)
print("Amount of regions of interest: " +  str(counter_output))
print("--------------- FINISHED REGIONOFINTERESTFILTER.PY ----------------")