library(BSgenome)
library(Biostrings)
library(seqinr)

chr_19 <- read.fasta(file = "chr19.fa.gz", seqtype = "DNA", forceDNAtolower = TRUE)
sub_len <- 5000000
sub_idx <- sample(1:(length(chr_19[[1]]) - sub_len + 1), 1)
sub <- chr_19[[1]][sub_idx:(sub_idx + sub_len - 1)]
write.fasta(sequences = sub, names = "partof19", file.out = "partof19.fa")

#SNVs
snv_indices <- sample(1:length(sub), 400)
snv_ref_nucs <- sub[snv_indices]
snv_alt_nucs <- sample(c("a", "c", "g", "t"), 400, replace = TRUE)
sub[snv_indices] <- snv_alt_nucs

bio_sub <- BString(paste(sub, collapse = ""))
#Insertion
ins_indices <- sample(0:(length(sub)-1), 300)
ins_ref_nucs <- vector()
ins_alt_nucs <- vector()
ins_positions <- vector()
for (ins_pos in ins_indices)
{
  if (ins_pos == 0)
  {
    ins_len <- sample(1:5, 1)
    ins_positions <- c(ins_positions, ins_len + 2)
    ins_ref_nucs <- c(ins_ref_nucs, toString(subseq(bio_sub, start = 1, end = 1)))
    subseq(bio_sub, start = 1, end = 1) <- BString(paste(c(toString(bio_sub[1]),sample(c("a", "c", "g", "t"), ins_len, replace = TRUE)), collapse = ""))
    ins_alt_nucs <- c(ins_alt_nucs, toString(subseq(bio_sub, start = 1, end = 1+ins_len)))
  }
  else
  {
    ins_len <- sample(1:5, 1)
    ins_positions <- c(ins_positions, ins_pos)
    ins_ref_nucs <- c(ins_ref_nucs, toString(subseq(bio_sub, start = ins_pos, end = ins_pos+1)))
    subseq(bio_sub, start = ins_pos+1, end = ins_pos+1) <- BString(paste(c(toString(bio_sub[ins_pos+1]), sample(c("a", "c", "g", "t"), ins_len, replace = TRUE)), collapse = ""))
    ins_alt_nucs <- c(ins_alt_nucs, toString(subseq(bio_sub, start = ins_pos, end = ins_pos + 1 + ins_len)))
  }
}

#Deletions
del_indices <- sample(0:(length(bio_sub)-1), 300)
del_ref_nucs <- vector()
del_alt_nucs <- vector()
del_positions <- vector()
for (del_pos in del_indices)
{
  if (del_pos >= 1 & del_pos <= (length(bio_sub)-5))
  {
    del_positions <- c(del_positions, del_pos)
    del_len <- sample(1:5, 1)
    del_ref_nucs <- c(del_ref_nucs, toString(subseq(bio_sub, start = del_pos, end = del_pos + del_len)))
    del_alt_nucs <- c(del_alt_nucs, toString(bio_sub[del_pos]))
    subseq(bio_sub, start = del_pos+1, end = del_pos+del_len) <- NULL
  }
  else if (del_pos >= (length(bio_sub)-4))
  {
    del_positions <- c(del_positions, del_pos)
    del_len <- sample(1:(length(bio_sub)-del_pos))
    del_ref_nucs <- c(del_ref_nucs, toString(subseq(bio_sub, start = del_pos, end = del_pos + del_len)))
    del_alt_nucs <- c(del_alt_nucs, toString(bio_sub[del_pos]))
    subseq(bio_sub, start = del_pos+1, end = del_pos+del_len) <- NULL
  }
  else if (del_pos == 0)
  {
    del_len <- sample(1:5, 1)
    del_positions <- c(del_positions, 1+del_len)
    del_ref_nucs <- c(del_ref_nucs, toString(subseq(bio_sub, start = 1, end = del_len+1)))
    del_alt_nucs <- c(del_alt_nucs, toString(bio_sub[del_len+1]))
    subseq(bio_sub, start = 1, end = del_len) <- NULL
  }
}

write.fasta(sequences = bio_sub, names = "sample", file.out = "Part1/Sample.fa")
CHROM <- rep(19, 1000)
POS <- c(snv_indices, ins_positions, del_positions)
REF <- c(snv_ref_nucs, ins_ref_nucs, del_ref_nucs)
ALT <- c(snv_alt_nucs, ins_alt_nucs, del_alt_nucs)
vcf <- data.frame(CHROM, POS, REF, ALT)
write.table(vcf, "Part1/Sample.vcf", sep = "\t", row.names = FALSE)
