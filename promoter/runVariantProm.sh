
# read in parameters
echo "name: $1"
echo "output_path: $2"
echo "input_path: $3"
echo "promoter_path: $4"
echo "start: $5"
echo "end: $6"

name=$1
output_path=$2
input_path=$3
promoter_path=$4
start=$5
end=$6

# pattern matching for dna_vcf files
pattern="*.vcf"
files=$(find "$input_path" -name "$pattern")
full_output_path="$output_path/vcf_promoter_regions"

# create dir, check for output dir
if [ ! -d "$output_path"]; then
    mkdir "$output_path"
fi

# create subfolder
if [ ! -d "$full_output_path"]; then
    mkdir -p "$full_output_path"
fi

for file in $files; do
    full_input_path="$input_path/$file"
    echo filename "$file"
    # run script
    python3 "variantPromoterRegion.py" -n "$name" -o "$full_output_path" -v "$file" -p "$promoter_path" -s "$start" -e "$end"
done
# python3 variantPromoterRegion.py -n "alpha" -o "informationAndData/vcfsPromoter/" -v "informationAndData/test_vcf_file" -p informationAndData/GRCh37_promoterChrPos_testCopy.bed -s "500" -o "100"
#