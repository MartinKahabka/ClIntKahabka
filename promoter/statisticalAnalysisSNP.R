# author: Martin Kahabka
# conducts a fishers test on all variants and saves the results
print("--- START PROGRAMM STATISTICALANALYSISPROMOTER.R")
tryCatch({
  library(tzdb, lib.loc = "./rLibs")
  library(readr, lib.loc = "./rLibs")
  library(vroom, lib.loc = "./rLibs")
}, error = function(error) {
  library(readr)
})

# read in args
args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
print(input_path)
print(output_path)

# read in data, get total number of patients
data <- read.table(input_path, sep = "\t", header = TRUE)
num_variants <- nrow(data)

# p-value bonferroni correction
standart_p_value <- 0.05
corrected_p_value <- standart_p_value / num_variants

# fishers excat test
# dataframe chr pos sig? p-value
col_result <- c("ChrNum", "PosOnChr", "resultAnalysis", "p_value")
results_matrix <- matrix(nrow = 0, ncol = length(col_result))
for (i in seq_len(nrow(data))) {
  current <- data[i, ]
  # define matrix for fishers test
  values <- data.frame(
    "severe" = c(current[[3]], current[[4]]),
    "notSevere" = c(current[[5]], current[[6]]),
    row.names = c("variant", "not variant")
  )
  # compute fisher exact test for variant
  test <- fisher.test(values)
  # get values for ouput file
  p_value_of_test <- test$p.value
  if (p_value_of_test <= corrected_p_value) {
    result <- "significant"
  } else {
    result <- "not_significant"
  }
  # saves to matrix
  new_row <- c(current[[1]], current[[2]], result, p_value_of_test)
  results_matrix <- rbind(results_matrix, new_row)
}
# write over to data frame
results_analysis <- data.frame(
  results_matrix,
  row.names = seq_len(nrow(data)),
  stringsAsFactors = FALSE
)
colnames(results_analysis) <- col_result

# save in file
write_tsv(results_analysis, path = output_path)

print("--- SUCCESFUL ---")