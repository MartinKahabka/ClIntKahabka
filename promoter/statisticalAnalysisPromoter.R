args <- commandArgs(trailingOnly = TRUE)
path <- args[1]

# read in data, get total number of patients
data <- read.table(path, sep = "\t", header = TRUE)
num_variants <- nrow(data)

# p-value bonferroni correction
standart_p_value <- 0.05
corrected_p_value <- standart_p_value / num_variants

# fishers excat test
# input in for of
# sample corresponds to severe/not severe
# contingency table
#            | severe  | not severe
# variant    |  row[3] |    row[4]
# no variant |  row[5] |    row[6]
data_frame <- data.frame(column1 = c(1, 2, 3), column2 = c("a", "b", "c"))
print(data_frame)

for (i in seq_len(nrow(data))) {
    current <- data[i, ]
    # define matrix for fishers test
    values <- data.frame(
        "severe" = c(current[[3]], current[[4]]),
        "notSevere" = c(current[[5]], current[[6]]),
        row.names = c("variant", "not variant")
    )
    test <- fisher.test(values)
}