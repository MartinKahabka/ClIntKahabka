#!/bin/bash -v 3.2  # Use version 3.2

echo "name: $1"
echo "output_path: $2"
echo "input_path_vcf_patient: $3"
echo "gene_names: $4"
echo "database: $5"
echo "start: $6"
echo "end: $7"
echo "path_patient_info: $8"
echo "fast: $9"

name=$1
output_path=$2
input_path_vcf_patient=$3
gene_names=$4
database=$5
start=$6
end=$7
path_patient=$8
fast=$9

# path of output folder
full_output_path="$output_path/$name"

# path from filter
output_path_filtered="$full_output_path/regionsOfInterest.bed"

# path of output in variant promoter
output_path_vcf="$full_output_path/vcf_promoter_regions"
output_sum_path="$full_output_path/variants_per_promoter"

# part of output in summary for result
output_path_summary="$full_output_path/summary_promoter.tsv"
output_path_sum="$full_output_path/summary_sum_variants.tsv"

# part of output in statistical analysis
output_path_statistcal_snp="$full_output_path/statistical_result.tsv"
output_path_statistcal_sum="$full_output_path/statistical_sum_variant_result.tsv"

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

# filter from database
python3 regionOfInterestFilter.py -d "$database" -g "$gene_names" -o "$output_path_filtered"

# run variant promoter
sh ./runVariantFilter.sh $name $full_output_path $input_path_vcf_patient $output_path_filtered $start $end $fast $output_sum_path

# summary for result 
python3 summaryPromoterResults.py -i "$output_path_vcf" -p "$path_patient" -o "$output_path_summary"
python3 summarySumOfVariants.py -i "$output_sum_path" -p "$path_patient" -o "$output_path_sum"


# statistical analysis
Rscript statisticalAnalysisPromoter.R "$output_path_summary" "$output_path_statistcal_snp"
Rscript statisticalAnalysisSum.R "$output_path_sum" "$output_path_statistcal_sum"

