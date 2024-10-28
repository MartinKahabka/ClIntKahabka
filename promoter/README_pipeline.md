 README for bash script pipeljneScript.sh (pipeline PiaVaRR). For a more detailed view of the input file structures, see below

**parameter:** Script has 9 parameters `<br/>`
name: unique identifier for a pipeline run `<br/>`
output_path: path where the output folder is created `<br/>`
input_path_vcf_patient: vcf files of the patients to be analysed `<br/>`
gene_names: names of the specific gene to be analysed `<br/>`
database: database from with the regions are filtered `<br/>`
start: defines range of search downstream from start of region `<br/>`
end: defines range of search upstream from start of region `<br/>`
path_patient_info: information about patients (severe/not severe covid) `<br/>`
binary: if true binarySearch algorithm will be executed, if not linear sorted lists (linear) algorithm will be executed `<br/>`

**example**
"sh pipelineScript.sh test_run "informationAndData/output_promoter" "informationAndData/input" "informationAndData/feature.txt" "data/silencer_blood_predicted_hg19.bed" 500 100 "informationAndData/Q001H_sample_preparations_20230803115337_adjusted.tsv" true none"

**pipeline PiaProVa**`<br/>`
The PiaProVa pipeline is created to analyse DNA mutations from given patients data (.vcf files). The pipeline takes a dataset containing gene regions (promoter/enhancer/silencer) and compares this data to the specified genes to filter for "regions of interest".
Each region in the datasets has a gene associated with it, if this gene is in the list of specified genes, the region will be consider a "region of interest". `<br/>`
**Currently there are two mode of statistical analysis:**`<br/>`
**SNP:** Indiviudual mutations are counted in patients with severe and not severe covid conditions and for each position with a mutation a fisher exact test is used to test for significant differences between patients with severe and not severe covid conditions.`<br/>`
**Quantity of variant:** The number of mutations for each patients per region of interest is counted and the collected. A KS-test is used to compare the sums from patients with severe covid to patients with not severe conditions.`<br/>`

**pipeline structure**`<br/>`
see file "structure.drawio" (draw.io can be used to view the file)

**runtime:**`<br/>`
The rumtime differs for both algorithms
**Sorted lists (linear):** O(n)
**binary search:** O(log(n) * (p * r) `<br/>`
where`<br/>`
n: number of variants of patient `<br/>`
p: number of regions of interest `<br/>`
r: range of regions of interest `<br/>`

**important notes:**

1. If the pipeline executes both statistical analyses the p-value correction is not done collectively, but in both analyses seperately. The used p-value correction is the bonferroni correction. This p-value correction is outdated. Use own p-value correction `<br/>`
2. If no region is consider a region of interest, meaning every region is filtered out and no region is analysed further, the pipeline can crash. `<br/>`
3. Only files that follow the pattern r"FO\d*x\d*" are read in `<br/>`

**Structure of input files**`<br/>`
If new data is added/other data is used, the pipeline only works correctly if files are structure as follows:`<br/>`
**input_path_vcf_patient:** vcf files of the patients to be analysed:`<br/>`
see formal structure of .vcf files (https://en.wikipedia.org/wiki/Variant_Call_Format)`<br/>`

**gene_names:** names of the specific gene to be analysed`<br/>`
.txt file with names of genes sepeated by a newline: `<br/>`
Example`<br/>`
"""`<br/>`
CD177`<br/>`
WFDC1`<br/>`
CD177P1`<br/>`
...`<br/>`
"""`<br/>`

**database:** database from with the regions are filtered `<br/>`
.bed files. See /promoter/data for example of structure

**path_patient_info:** information about patients (severe/not severe covid) `<br/>`
.tsv file where every row contains information about one patient. The patients .vcf is identified by an id in the second column.`<br/>`
Example: .vcf file: FO13636x01.vcf  | id: FO13636x01

The codition of the patients are read in from column 16.

For more details see summarySumOfVariants.py/summaryPromoterResults.py, function add_conditions
