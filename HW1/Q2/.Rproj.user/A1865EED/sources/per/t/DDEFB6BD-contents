library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(genemodel)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#A
#symbol = readline(prompt = "Enter the gene symbol: ")
symbol <- "RPL5"
gene_id <- select(org.Hs.eg.db, keys = symbol, keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID
txs <-select(txdb, keys = gene_id, keytype = "GENEID", columns = c("TXID", "TXNAME"))
print(txs)
#B
chr <- select(txdb, keys = gene_id, keytype = "GENEID", columns = "EXONCHROM")$EXONCHROM
chr_length <- seqlengths(Hsapiens)[[chr]]
exons <- select(txdb, keys = gene_id, keytype = "GENEID", columns = c("EXONSTART", "EXONEND"))
starts <- exons[["EXONSTART"]]
ends <- exons[["EXONEND"]]
start <- min(starts)
end <- max(ends)
gene_length <- end - start + 1
print(gene_length/chr_length)
#########################################
#Genemodel Plot
tx_names <- txs$TXNAME
exons_all <- exonsBy(txdb, "tx", use.names = TRUE)[tx_names]
introns_all <- intronsByTranscript(txdb, use.names = TRUE)[tx_names]
for (i in length(exons_all):1)
{
  exons <- exons_all[i]
  introns <- introns_all[i]
  type <- vector()
  exon_range <- ranges(exons)[[1]]
  utr_5 <- start(ranges(exons)[[1]][1])
  bpstart <- utr_5 - 100
  start(exon_range) <- start(exon_range) - bpstart
  end(exon_range) <- end(exon_range) - bpstart
  exons_interval <- as.character(exon_range)
  type <- rep("coding_region", length(exons_interval))
  intron_range <- ranges(introns)[[1]]
  start(intron_range) <- start(intron_range) - bpstart
  end(intron_range) <- end(intron_range) - bpstart
  intron_interval <- as.character(intron_range)
  type <- c(type, rep("intron", length(intron_interval)))
  utr_5 <- "1-100"
  type <- c(type, "5' utr")
  utr_3 <- end(ranges(exons)[[1]][length(ranges(exons)[[1]])]) - bpstart
  utr_3 <- paste(c(utr_3, utr_3+10), collapse = "-")
  type <- c(type, "3' utr")
  coordinates <- c(exons_interval, intron_interval, utr_5, utr_3)
  data <- data.frame(type, coordinates)
  bpstop <- end(ranges(exons)[[1]][length(ranges(exons)[[1]])]) + 100
  genemodel.plot(model = data, start = bpstart, bpstop = bpstop, orientation = "reverse")
}
