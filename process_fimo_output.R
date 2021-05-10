library(data.table)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)
fasta_file <- args[1]
fimo_tsv <- args[2]
out_tsv <- args[3]

test <- readDNAStringSet(fasta_file,"fasta")
a <- fread(fimo_tsv, header=T)
motifs <- unique(a$motif_id)
temp <- split(a, f=factor(a$motif_id))

# p-value matrix
motif_p_m <- matrix(1, length(test), length(motifs), dimnames=list(names(test), motifs))
for (i in 1:length(temp)){
  motif <- names(temp)[i]
  pvalue_tmp <- tapply(temp[[i]]$'p-value', temp[[i]]$sequence_name, min) # for each peak, find the min pvalue
  motif_p_m[names(pvalue_tmp), motif] <- pvalue_tmp
  print(i)
}

fwrite(motif_p_m, file=out_tsv)

