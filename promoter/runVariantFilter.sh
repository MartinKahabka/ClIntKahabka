#!/bin/bash -v 3.2  # Use version 3.2

echo "--------------- START RUNVARIANTPROM.SH ----------------"
echo "Input arguments"
# read in parameters
echo "name: $1"
echo "output_path: $2"
echo "input_path: $3"
echo "region_of_interest_path: $4"
echo "start: $5"
echo "end: $6"
echo "binary_algorithm: $7"
echo "output_regions_path: $8"

name=$1
output_path=$2
input_path=$3
promoter_path=$4
start=$5
end=$6
fast=$7
output_sum_path=$8

# pattern matching for dna_vcf files
pattern="*.vcf"
files=$(find "$input_path" -name "$pattern")
full_output_path="$output_path/snp_in_regions"

# create dir, check for output dir
if [ ! -d "$output_path" ]; then
    mkdir "$output_path"
fi

# create subfolder
if [ ! -d "$full_output_path" ]; then
    mkdir -p "$full_output_path"
fi

if $fast == "true"; then
    for file in $files; do
    full_input_path="$input_path/$file"
    echo filename: "$file"
    # run script
    python3 "variantFilter_binary.py" -n "$name" -o "$full_output_path" -v "$file" -p "$promoter_path" -s "$start" -e "$end" -f "$output_sum_path"
    done
else
    for file in $files; do
        full_input_path="$input_path/$file"
        echo filename "$file"
        # run script
        python3 "variantFilter_linear.py" -n "$name" -o "$full_output_path" -v "$file" -p "$promoter_path" -s "$start" -e "$end" -f "$output_sum_path"
    done
fi
# python3 variantPromoterRegion.py -n "alpha" -o "informationAndData/vcfsPromoter/" -v "informationAndData/test_vcf_file" -p informationAndData/GRCh37_promoterChrPos_testCopy.bed -s "500" -e "100"

echo "--------------- FINISHED RUNVARIANTPROM.SH ----------------"