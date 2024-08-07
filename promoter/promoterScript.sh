# folder structure
# output_promoter/name1/vcf_promoter_regions/...
# output_promoter/name1/summary
# output_promoter/name1/statistical
#               ./name2
#!/bin/bash -v 3.2  # Use version 3.2

echo "name: $1"
echo "output_path: $2"
echo "input_path_vcf_patient: $3"
echo "promoter_path: $4"
echo "start: $5"
echo "end: $6"
echo "path_patient_info: $7"
echo "fast: $8"
echo "mode: $9"

name=$1
output_path=$2
input_path_vcf_patient=$3
promoter_path=$4
start=$5
end=$6
path_patient=$7
fast=$8
mode=$9

# path of output folder
full_output_path="$output_path/$name"

# path of output in first part
output_path_vcf="$full_output_path/vcf_promoter_regions"
output_sum_path="$full_output_path/variants_per_promoter"

# part of output in second part
output_path_summary="$full_output_path/summary_promoter.tsv"
output_path_sum="$full_output_path/summary_sum_variants.tsv"

# part of output in third part
output_path_statistcal="$full_output_path/statistical_result.tsv"
output_path_statistcal="$full_output_path/statistical_sum_variant_result.tsv"

# check and create necessary output folder
if [ ! -d "$output_path" ]; then
    mkdir "$output_path"
fi

if [ ! -d "$full_output_path" ]; then
    mkdir -p "$full_output_path"
fi

if [ ! -d "$output_path_vcf" ]; then
    mkdir -p "$output_path_vcf"
fi

if [ ! -d "$output_sum_path" ]; then
    mkdir -p "$output_sum_path"
fi

# run first part
sh ./runVariantProm.sh $name $full_output_path $input_path_vcf_patient $promoter_path $start $end $fast $output_sum_path

# run second part
python3 summaryPromoterResults.py -i "$output_path_vcf" -p "$path_patient" -o "$output_path_summary"
python3 summarySumOfVariants.py -i "$output_sum_path" -p "$path_patient" -o "$output_path_sum"


# run third part
Rscript statisticalAnalysisPromoter.R "$output_path_summary" "$output_path_statistcal"
Rscript statisticalAnalysisSum.R "$output_path_sum" "$output_path_statistcal"

