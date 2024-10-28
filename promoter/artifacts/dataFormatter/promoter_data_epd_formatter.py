# author martin kahabka
# example input line: 
# chr1	894625	894685	NOC2L_1	900	-	894625	894636

# example output line:
# chr1	894625	894685	NOC2L_1 1 -
with open("artifacts/promoter_data_EPDnew_hg19.bed", "r") as epd_file, open("data/promoter_data_EPDnew_hg19_format.bed", "w") as formatted_file:
    for promoter in epd_file:
        # get contents of line, \t seperated values
        promoter = promoter.strip()
        content = promoter.split("\t")
        
        chr_num = content[0]
        start = int(content[1])
        end = int(content[2])
        # name of gene
        name = content[3]
        value = 1
        strand_dir = content[5]
        
        # combine to string + newline
        data_string = "\t".join([chr_num, str(start), str(end), name, str(value), strand_dir]) +  "\n"
        
        # write to file
        formatted_file.write(data_string)
        

        
        