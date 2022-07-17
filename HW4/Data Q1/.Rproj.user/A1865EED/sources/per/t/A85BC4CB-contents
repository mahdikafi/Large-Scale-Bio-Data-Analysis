library(edgeR)

count.data <- read.table(file = 'GSE60450_Lactation-GenewiseCounts.txt', header = TRUE, sep = "\t")
n_col <- dim(count.data)[2]
# CPM Normalization
cpm.normalized <- cpm(y = count.data[, 3:n_col], log = TRUE)
write.csv(cpm.normalized, file = "cpm_normalized.csv", row.names = count.data[, 1])

# TPM Normalization 
count.data <- read.table(file = 'GSE60450_Lactation-GenewiseCounts.txt', header = TRUE, sep = "\t")
count.data[, 3:n_col] <- count.data[, 3:n_col] / (count.data[, 2] / 1000)
count.data[, 3:n_col] <- count.data[, 3:n_col] / (colSums(count.data[, 3:n_col])/1000000)
count.data[, 3:n_col] <- log2(count.data[, 3:n_col]+0.25)
write.csv(count.data[, 3:n_col], file = "tpm_normalized.csv", row.names = FALSE)
 
#DESeq2 Normalization has been done using Python. 


