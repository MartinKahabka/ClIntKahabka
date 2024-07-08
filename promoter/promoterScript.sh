# folder structure
# output_promoter/name1/vcf_promoter_regions/...
# output_promoter/name1/summary
# output_promoter/name1/statistical
#               ./name2

echo "name: $1"
echo "output_path: $2"
echo "input_path_vcf_patient: $3"
echo "promoter_path: $4"
echo "start: $5"
echo "end: $6"

name=$1
output_path=$2
input_path_vcf_patient=$3
promoter_path=$4
start=$5
end=$6

full_output_path="$output_path/$name"
# create output_path
if [ ! -d "$output_path"]; then
    mkdir "$output_path"
fi

if [ ! -d "$full_output_path"]; then
    mkdir "$full_output_path"
fi