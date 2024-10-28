# author: Martin Kahabka
# conducts a fishers test on all variants and saves the results
print("--------------- START STATISTICALANALYSISSNP.R ---------------")
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
print(paste("Input path:", input_path))
print(paste("Output path:", output_path))

# read in data, get total number of patients
data <- read.table(input_path, sep = "\t", header = TRUE)
num_variants <- nrow(data)

# p-value
standart_p_value <- 0.05

# fishers excat test
# dataframe chr pos sig? p-value
col_result <- c("ChrNum", "PosOnChr", "p_value")
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
  # saves to matrix
  new_row <- c(current[[1]], current[[2]], test$p.value)
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
write_tsv(results_analysis, file = output_path)

print("--------------- FINISHED STATISTICALANALYSISSNP.R ---------------")