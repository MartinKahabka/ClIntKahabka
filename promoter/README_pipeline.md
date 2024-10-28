 README for bash script pipeljneScript.sh (pipeline PiaVaRR). For a more detailed view of the input file structures, see below

**parameter:** Script has 9 parameters  
name: unique identifier for a pipeline run  
output_path: path where the output folder is created  
input_path_vcf_patient: vcf files of the patients to be analysed  
gene_names: names of the specific gene to be analysed  
database: database from with the regions are filtered  
start: defines range of search downstream from start of region  
end: defines range of search upstream from start of region  
path_patient_info: information about patients (severe/not severe covid)  
binary: if true binarySearch algorithm will be executed, if not linear sorted lists (linear) algorithm will be executed  

**example**  
"sh pipelineScript.sh test_run "informationAndData/output_promoter" "informationAndData/input" "informationAndData/feature.txt" "data/silencer_blood_predicted_hg19.bed" 500 100 "informationAndData/Q001H_sample_preparations_20230803115337_adjusted.tsv" true none"

**pipeline PiaProVa**  
The PiaProVa pipeline is created to analyse DNA mutations from given patients data (.vcf files). The pipeline takes a dataset containing gene regions (promoter/enhancer/silencer) and compares this data to the specified genes to filter for "regions of interest".  
Each region in the datasets has a gene associated with it, if this gene is in the list of specified genes, the region will be consider a "region of interest".  
**Currently there are two mode of statistical analysis:**  
**SNP:** Indiviudual mutations are counted in patients with severe and not severe covid conditions and for each position with a mutation a fisher exact test is used to test for significant differences between patients with severe and not severe covid conditions.  
**Quantity of variant:** The number of mutations for each patients per region of interest is counted and the collected. A KS-test is used to compare the sums from patients with severe covid to patients with not severe conditions.  

**pipeline structure**  
see file "structure.drawio" (draw.io can be used to view the file)

**runtime:**  
The rumtime differs for both algorithms  
**Sorted lists (linear):** O(n)  
**binary search:** O(log(n) * (p * r)  
where  
n: number of variants of patient  
p: number of regions of interest  
r: range of regions of interest  

**important notes:**

1. If the pipeline executes both statistical analyses the p-value correction is not done collectively, but in both analyses seperately. The used p-value correction is the bonferroni correction. This p-value correction is outdated. Use own p-value correction.  
2. If no region is consider a region of interest, meaning every region is filtered out and no region is analysed further, the pipeline can crash.  
3. Only files that follow the pattern r"FO\d*x\d*" are read in.  

**Structure of input files**`<br/>`
If new data is added/other data is used, the pipeline only works correctly if files are structure as follows:  
**input_path_vcf_patient:** vcf files of the patients to be analysed:  
see formal structure of .vcf files (https://en.wikipedia.org/wiki/Variant_Call_Format)  

**gene_names:** names of the specific gene to be analysed  
.txt file with names of genes sepeated by a newline:  
Example  
"""  
CD177  
WFDC1  
CD177P1  
...  
"""  

**database:** database from with the regions are filtered  
.bed files. See /promoter/data for example of structure  

**path_patient_info:** information about patients (severe/not severe covid)  
.tsv file where every row contains information about one patient. The patients .vcf is identified by an id in the second column.  
Example: .vcf file: FO13636x01.vcf  | id: FO13636x01  

The codition of the patients are read in from column 16.  

For more details see summarySumOfVariants.py/summaryPromoterResults.py, function add_conditions  
