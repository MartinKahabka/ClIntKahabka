import re
# author martin kahabka

# example output line:
# chr1	894625	894685	NOC2L_1 1 -
names = {}
with open("artifacts/enhancer_pbmc_interaction_hg19.txt", "r") as input_file, open("data/enhancer_pbmc_interaction_hg19.bed", "w") as output_file:
    for line in input_file:
        line = line.strip()
        
        # splits the line at the given delimiter
        information_region =  re.split(':|_|\$|\t', line)
        information_range = re.split('-', information_region[1])
        print(information_region)
        print(information_range)
        chr_num = information_region[0]
        start = int(information_range[0])
        end = int(information_range[1])
        value = 1
        strand_dir = information_region[6]
        print(end - start)
        # name of gene, convert to name of region (name -> name_x)
        name_gene = information_region[3]
        if name_gene in names:
            names[name_gene] += 1
        else:
            names[name_gene] = 1
        # name of region
        name = name_gene + "_" + str(names[name_gene])
    
        # combine to string + newline
        data_string = "\t".join([chr_num, str(start), str(end), name, str(value), strand_dir]) +  "\n"
        
        # write to file
        # output_file.write(data_string)