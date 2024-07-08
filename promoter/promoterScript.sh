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
echo "path_patient_info: $7"

name=$1
output_path=$2
input_path_vcf_patient=$3
promoter_path=$4
start=$5
end=$6
path_patient=$7

full_output_path="$output_path/$name"
output_path_vcf="$full_output_path/vcf_promoter_regions"
output_path_summary="$full_output_path/summary_promoter.tsv"
output_path_statistcal="$full_output_path/statistical_result.tsv"
# create output_path
if [ ! -d "$output_path" ]; then
    mkdir "$output_path"
fi

if [ ! -d "$full_output_path" ]; then
    mkdir -p "$full_output_path"
fi

if [ ! -d "$output_path_vcf" ]; then
    mkdir -p "$output_path_vcf"
fi

# run first script
sh ./runVariantProm.sh $name $full_output_path $input_path_vcf_patient $promoter_path $start $end

# run second script
python3 summaryPromoterResults.py -i "$output_path_vcf" -p "$path_patient" -o "$output_path_summary"

# run third script
Rscript statisticalAnalysisPromoter.R "$output_path_summary" "$output_path_statistcal"