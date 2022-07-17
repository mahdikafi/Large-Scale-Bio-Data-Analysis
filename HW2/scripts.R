library("BSgenome")
library("Biostrings")
library("GenomicRanges")

#PART0: basics 
#GOAL1: write a function to generate a random genome with a given length
ref_length <- 1000000
generate_genome <- function(genome_length)
{
  genome <- DNAString(paste0(sample(c("A", "C", "G", "T"), size = genome_length, replace = TRUE), collapse = ""))  
  return(genome)
}
ref <- generate_genome(ref_length)
print(ref)
#GOAL2: write a function to generate random reads with given length and numbers
read_length <- 150
read_num <- 10000
generate_reads <- function(ref, read_length, read_num)
{
  starts <- sample(1:(ref_length - read_length + 1), size = read_num, replace = TRUE)
  REF <- DNAStringSet(ref)
  names(REF) <- "ref"
  gr <- GRanges(seqnames = "ref", ranges = IRanges(start = starts, width = read_length), strand = rep("+", read_num))
  reads <- getSeq(REF, gr)
  names(reads) <- paste0("Read", 1:read_num)
  return(reads)  
}
naive_reads <- generate_reads(ref, read_length, read_num)
print(naive_reads)
#GOAL3: write simulated genome and reads in a file
REF <- DNAStringSet(ref)
names(REF) <- "ref"
writeXStringSet(x=REF, filepath = "ref.fasta")
writeXStringSet(x=naive_reads, filepath = "reads.fastq", format = "fastq")
#PART1: simulate a dna fragmentation with sonication 
#Write a function in R with the following inputs:
#1- Genome
#2- Number of Fragments
#3- Mean of Fragments
#4- Variance of Fragments
#The output is a list of fragments sampled from a normal distribution with specified strands
#GOAL1: Plot a histogram of fragment lengths. 
simulate_sonicator <- function(ref, fragment_num, fragment_width_mean, fragment_width_var)
{
  fragment_width <- floor(rnorm(fragment_num, fragment_width_mean, sqrt(fragment_width_var)))
  fragment_starts <- sample(1:(ref_length - max(fragment_width) + 1), size = fragment_num, replace = FALSE)
  REF <- DNAStringSet(ref)
  names(REF) <- "ref"
  gr <- GRanges(seqnames = "ref", ranges = IRanges(start = fragment_starts, width = fragment_width), strand = sample(c("+", "-"), size = fragment_num, replace = TRUE))
  fragments <- getSeq(REF, gr)
  return(fragments)
}
fragment_num <- read_num*10
fragment_width_mean <- 500
fragment_width_var <- 40
fragments <- simulate_sonicator(ref, fragment_num, fragment_width_mean, fragment_width_var)
hist(width(fragments), breaks = 100)
print(fragments)

#PART2: simulate PCR
#Write a function in R with the following inputs:
#1- output of fragments
#2- Lambda_min
#3- Lambda_max
#4- modes
#The output is the number of copies for each fragment.
#Modes accepts two cases: in case one the copies are drawn from a poisson distribution with lambda_min.
#In case two the copies are sampled from poisson distribution with lambdas uniformly selected from Lambda_min and Lambda_max.
#GOAL1: Predict the complexity of a given library from N samples from single end reads
#GOAL2: Predict the complexity of a given library from N samples from paired end reads
PCR <- function(fragments, lambda_min, lambda_max=NULL, mode=1){
  if (mode == 1){
    return(rpois(length(fragments), lambda = lambda_min))
  }
  if (mode == 2){
    return(rpois(length(fragments), lambda = sample(x = lambda_min:lambda_max, size = length(fragments), replace = TRUE)))
  }
}

cnts <- PCR(fragments, lambda_min = 10, mode = 1)

single_end_sequencing <- function(fragments, cnts, N, L){
  library <- rep(fragments, cnts)
  lib.indices <- rep(1:length(fragments), cnts)
  read.indices<- sample(1:length(library), size = N)
  selected_fragments <- library[read.indices]
  names(selected_fragments) <- paste0("s_fragments", 1:N)
  gr <- GRanges(seqnames = paste0("s_fragments", 1:N), ranges = IRanges(start = 1, width = L), strand = rep("+", N))
  reads <- getSeq(selected_fragments, gr)
  indices <- lib.indices[read.indices]
  metadata(reads)$indices <- indices
  return(reads)
}

reads <- single_end_sequencing(fragments, cnts, read_num, read_length)
indices <- metadata(reads)$indices
M <- length(unique(indices))

C <- seq(M, 10*read_num, 1000)
M_hat <- C*(1-exp(-read_num/C))
plot(C, M_hat, 'l')

estimated_complexity <- C[which.min((M_hat-M)^2)]
print(paste0("estimated complexity: ",estimated_complexity))
print(paste0("number of fragments: ", fragment_num))
print((fragment_num-estimated_complexity)/fragment_num)

#GOAL3: Plot Coverage vs number of samples from single end reads
no_reads <- seq(0, read_num*100, 1000)
independent_fragments <- estimated_complexity*(1-exp(-no_reads/estimated_complexity))
plot(no_reads, independent_fragments, 'l')

#GOAL4: Predict if adding more reads with increase the coverage
max_coverage <- max(independent_fragments)
print(fragment_num - max_coverage)
