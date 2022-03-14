library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
source("~/myscripts/NGS_scripts/get_seqs.R")

gr <- import.bed('gr_1344.bed')
seqs <- get.seqs(BSgenome.Hsapiens.UCSC.hg19, gr)
writeXStringSet(seqs, "gr.fasta", format="fasta")

system("fimo /home/yuanh/programs/source/motif_databases.12.18/CIS-BP/Homo_sapiens.meme gr.fasta")
