library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Ecoli.NCBI.20080805)
library(RSVSim)

#PART0: basics 
#We work on three genomes all having same size as ECOLI
#GENOME 1: a random genome of  ecoli genome size
#GENOME 2: ecoli
#GENOME 3: a part of human chromosome 20 starting from 20Mbp and having the same size as ecoli genome.
#3 - Write a program that assembles genomes from a given set of reads ( you may only call available packages such as velvet).
#1-  write a program and prepare all these genomes
ecoli <- Ecoli$NC_000913
ECOLI <- DNAStringSet(ecoli)
names(ECOLI) <- "Ecoli"
writeXStringSet(ECOLI, filepath = "data/ecoli.fasta")
print(ECOLI)
len <- length(ecoli)
ran_genome <- DNAString(paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = ""))
RAN_GENOME <- DNAStringSet(ran_genome)
names(RAN_GENOME) <- "RandomGenome"
writeXStringSet(RAN_GENOME, filepath = "data/random_genome.fasta")
print(RAN_GENOME)
human <- Views(subject = Hsapiens$chr20, start = 20000000, width = len)[[1]]
HUMAN <- DNAStringSet(human)
names(HUMAN) <- "Human"
writeXStringSet(HUMAN, filepath = "data/human.fasta")
print(HUMAN)
#2- Write a program that generates reads from a genome with a given error rate.
#(Consider only substitution errors)
read_length <- 150
read_num <- 10000
generate_reads <- function(ref, read_length, read_num, error_rate, read_name)
{
  starts <- sample(1:(width(ref) - read_length + 1), size = read_num, replace = TRUE)
  gr <- GRanges(seqnames = names(ref), ranges = IRanges(start = starts, width = read_length), strand = sample(c("+", "-"), size = read_num, replace = TRUE))
  reads <- getSeq(ref, gr)
  print(length(reads))
  # reads <- Views(ref, start = starts, width = read_length)
  indices <- rep(FALSE, length(reads) * width(reads)[1] )
  #error_rate <- 1/1000
  indices[sample(1:(length(reads) * width(reads)[1]), size = ceiling(length(reads) * width(reads)[1] * error_rate))] <- TRUE
  at <- matrix(data = indices, nrow = length(reads), ncol = width(reads)[1])
  letter_subject <- DNAString(paste0(sample(c("A", "C", "G", "T"), size = ceiling(length(reads)*width(reads)[1]*error_rate), replace = TRUE), collapse = ""))
  starts <- sample(1:(ceiling(length(reads)*width(reads)[1]*error_rate) - max(rowSums(at))), size = length(reads), replace = TRUE)
  letter <- as(Views(letter_subject, start = starts, width = rowSums(at)), "DNAStringSet")
  subs_reads <- replaceLetterAt(reads, at, letter)
  names(subs_reads) <- paste0("Read", 1:read_num)
  writeXStringSet(x = subs_reads, filepath = paste0(c("data/", read_name, ".fastq"), collapse = ""), format = "fastq")
  return(subs_reads)  
}
error_rate <- 1/1000
ecoli_reads <- generate_reads(ECOLI, read_length, read_num, error_rate, "ecoli_reads")
ecoli_reads

human_reads <- generate_reads(HUMAN, read_length, read_num, error_rate, "human_reads")
human_reads

random_reads <- generate_reads(RAN_GENOME, read_length, read_num, error_rate, "random_reads")
random_reads

# 3 - Write a program that assembles genomes from a given set of reads ( you may only call available packages such as velvet).
system("velveth ecoli_out/ 21 -fastq -short data/ecoli_reads.fastq")
system("velvetg ecoli_out/ -cov_cutoff auto")
system("velveth human_out/ 21 -fastq -short data/human_reads.fastq")
system("velvetg human_out/ -cov_cutoff auto")
system("velveth random_out/ 21 -fastq -short data/random_reads.fastq")
system("velvetg random_out/ -cov_cutoff auto")

# PART1: DENOVO Assembly 
# Compare the performance of the denovo assembly based on the following variables.
# 1- Genome size
# 2- Number of reads
# 3- Read length
# 4- Error rate in reads
# GOAL1: Plot N50 vs coverage depth.
n50 <- function(velvetg_output){
  str <- strsplit(strsplit(output[length(output)], ",")[[1]][1], " ")[[1]]
  n50 <- str[length(str)]
  return(strtoi(n50))
}

n50.vec <- c()
for (coverage in 1:20){
  read_num <- floor(coverage*len/read_length)
  reads <- generate_reads(ECOLI, read_length, read_num, error_rate, "ecoli_reads")
  system("velveth ecoli_out/ 21 -fastq -short data/ecoli_reads.fastq")
  output <- system("velvetg ecoli_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec <- c(n50.vec, n50.value)
}
plot(1:20, n50.vec, xlab= "coverage", ylab="n50", type="b", main="n50 vs coverage: ECOLI")

n50.vec.human <- c()
for (coverage in 1:20){
  read_num <- floor(coverage*len/read_length)
  reads <- generate_reads(HUMAN, read_length, read_num, error_rate, "human_reads")
  system("velveth human_out/ 21 -fastq -short data/human_reads.fastq")
  output <- system("velvetg human_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec.human <- c(n50.vec.human, n50.value)
}
plot(1:20, n50.vec.human, xlab= "coverage", ylab="n50", type="b", main="n50 vs coverage:HUMAN")

n50.vec.random <- c()
for (coverage in 1:20){
  read_num <- floor(coverage*len/read_length)
  reads <- generate_reads(RAN_GENOME, read_length, read_num, error_rate, "random_reads")
  system("velveth random_out/ 21 -fastq -short data/random_reads.fastq")
  output <- system("velvetg random_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec.random <- c(n50.vec.random, n50.value)
}
plot(1:20, n50.vec.random, xlab= "coverage", ylab="n50", type="b", main="n50 vs coverage:RANDOM")
# GOAL2: Plot N50 vs error rate. 
read_length <- 150
read_num <- 25000
n50.vec.ecoli.er <- c()
for (error_rate in seq(0.0001, 0.5, 0.005)){
  reads <- generate_reads(ECOLI, read_length, read_num, error_rate, "ecoli_reads")
  system("velveth ecoli_out/ 21 -fastq -short data/ecoli_reads.fastq")
  output <- system("velvetg ecoli_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec.ecoli.er <- c(n50.vec.ecoli.er, n50.value)
}
plot(seq(0.0001, 0.5, 0.005), n50.vec.ecoli.er, xlab= "error rate", ylab="n50", type="b", main="n50 vs error rate: ECOLI")

n50.vec.human.er <- c()
for (error_rate in seq(0.0001, 0.5, 0.005)){
  reads <- generate_reads(HUMAN, read_length, read_num, error_rate, "human_reads")
  system("velveth human_out/ 21 -fastq -short data/human_reads.fastq")
  output <- system("velvetg human_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec.human.er <- c(n50.vec.human.er, n50.value)
}
plot(seq(0.0001, 0.5, 0.005), n50.vec.human.er, xlab= "error rate", ylab="n50", type="b", main="n50 vs error rate: HUMAN")

n50.vec.random.er <- c()
for (error_rate in seq(0.0001, 0.5, 0.005)){
  reads <- generate_reads(RAN_GENOME, read_length, read_num, error_rate, "random_reads")
  system("velveth random_out/ 21 -fastq -short data/random_reads.fastq")
  output <- system("velvetg random_out/ -cov_cutoff auto", intern = TRUE)
  n50.value <- n50(output[length(output)])
  n50.vec.random.er <- c(n50.vec.random.er, n50.value)
}
plot(seq(0.0001, 0.5, 0.005), n50.vec.random.er, xlab= "error rate", ylab="n50", type="b", main="n50 vs error rate: RANDOM")
# GOAL3: compare the three genomes.

# PART2: RESEQUENCING
# Consider the three genomes as reference genomes and create three new genomes with some mutations inserted.
# Mutations include substitutions and short indels.
# Write a program that calls variations from reads and a reference genome. ( you may use short read mapping such as bowtie)
mutation <- function(genome, mutation.rate){
  #mutation.rate <- 1/1000
  mutation.size <- floor(mutation.rate*width(genome))
  no.mutation <- floor(mutation.size/10)
  sim <- simulateSV(genome = genome, dels = no.mutation/2, sizeDels = 10)
  random.genome.mutation <- paste0(sample(c("A", "C", "G", "T"), size = width(genome), replace = TRUE), collapse = "")
  sim <- DNAStringSet(c(as.character(sim[[1]]), random.genome.mutation))
  names(sim) <- c("sim", "random")
  sim <- simulateSV(genome = sim, ins = no.mutation/2, sizeIns = 10)
  sim <- sim[[1]]
  no.snp <- mutation.size
  at <- sample(1:length(sim), size = no.snp)
  snp <- sample(c("A", "C", "G", "T"), size = no.snp, replace = TRUE)
  sim <- replaceLetterAt(sim, at, snp)
  return(sim)
}
mutation.rate <- 1/1000
mutated.ecoli <- mutation(ECOLI, mutation.rate)
mutated.random <- mutation(RAN_GENOME, mutation.rate)
mutated.human <- mutation(HUMAN, mutation.rate)
M.ECOLI <- DNAStringSet(mutated.ecoli)
M.HUMAN <- DNAStringSet(mutated.human)
M.RANDOM <- DNAStringSet(mutated.random)
writeXStringSet(M.ECOLI, "mutated/m_ecoli.fasta")
writeXStringSet(M.HUMAN, "mutated/m_human.fasta")
writeXStringSet(M.RANDOM, "mutated/m_random.fasta")
