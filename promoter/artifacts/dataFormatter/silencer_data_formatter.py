import csv
# author martin kahabka

# example output line:
# chr1	894625	894685	NOC2L_1 1 -
names = {}
input_format = "csv"
input_name = "silencer_blood_validated_hg19"
output_format = "bed"
with open("artifacts/" + input_name + "." + input_format, "r") as input_file, open("data/" + input_name + "." + output_format, "w") as output_file:
    csv_input_file = csv.reader(input_file)
    for line in csv_input_file:
        if line[0] != "id":
            chr_num = line[2]
            start = int(line[3])
            end = int(line[4])
            value = 1
            strand_dir = line[14]
            
            # name of gene, convert to name of region (name -> name_x)
            name_gene = line[11]
            if name_gene in names:
                names[name_gene] += 1
            else:
                names[name_gene] = 1
            # name of region
            name = name_gene + "_" + str(names[name_gene])
            
            
            
            # combine to string + newline
            data_string = "\t".join([chr_num, str(start), str(end), name, str(value), strand_dir]) +  "\n"
            
            # write to file
            output_file.write(data_string)