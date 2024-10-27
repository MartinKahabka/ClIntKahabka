# author: Martin Kahabka
# conducts a KS-test for every region of interest, saves the results in file
print("--------------- START STATISTICALANALYSISROI.R ---------------")
tryCatch({
  library(tzdb, lib.loc = "./rLibs")
  library(readr, lib.loc = "./rLibs")
  library(vroom, lib.loc = "./rLibs")
}, error = function(error) {
  library(readr)
})
# read in parameters
args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
print(paste("Input path:", input_path))
print(paste("Output path:", output_path))

# get number of promoter
lines <- readLines(input_path)
num_promoter <- length(lines) / 3

# read in data
data <- file(input_path, "r")

# create output file and write command
output_file <- file(output_path)
c1 <- "# statistical analysis of sum of variants per promoter"
c2 <- "# test used: KS - Test"
writeLines(c(c1, c2), output_file)

# input of cols: chr1	1010	TEST	sum_pos sum_neg p_value significant?
col_result <- c("ChrNum","PosOnChr", "PromoterName", "sum_pos", "sum_neg", "p_value")

# create necessary vars
results_matrix <- matrix(nrow = 0, ncol = length(col_result))
line <- " "
while (length(line) > 0) {
  # read three lines from file, representing one promoter
  line <- readLines(data, n = 3)
  header <- strsplit(line[1], "\t")[[1]]
  # convert to integer
  pos <- lapply(strsplit(line[2], "\t"), as.integer)[[1]]
  neg <- lapply(strsplit(line[3], "\t")[1], as.integer)[[1]]
  # check for NA values at end of file and break
  if (length(neg) <= 1) {
    break
  } else {
    # perform test and add to list of tests
    t <- ks.test(pos, y = neg, alternative = "two.sided")
    # save result
    new_row <- c(header, t$p.value)
    results_matrix <- rbind(results_matrix, new_row)
  }
}
# write over to data frame
results_analysis <- data.frame(
  results_matrix,
  row.names = seq_len(num_promoter),
  stringsAsFactors = FALSE
)
colnames(results_analysis) <- col_result

# write out output file
write_tsv(results_analysis, file = output_path)

close(output_file)

print("--------------- FINISHED STATISTICALANALYSISROI.R ---------------")