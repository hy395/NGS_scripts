# process IDR data
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
genome <- args[2]
p.value <- as.numeric(args[3])

a <- read.table(input, sep="\t", header=F, stringsAsFactors=F, colClasses="character")
gr <- GRanges(seqnames=a$V1, ranges=IRanges(as.numeric(a$V2), as.numeric(a$V3)), strand="*", 
              score=as.numeric(a$V5), signalValue=as.numeric(a$V7), summit=as.numeric(a$V10), localIDR=as.numeric(a$V11), globalIDR=as.numeric(a$V12))
thres <- -125*log2(p.value)

gr <- gr[seqnames(gr)%in%paste0("chr",c(1:22,"X","Y"))]
gr <- gr[gr$score>thres]

blacklist <- NULL
if (genome == "hg19") {
 blacklist <- import("/home/hy395/programs/genomes/blacklists/hg19.blacklist.bed.gz")
}
if (genome == "hg38") {
 blacklist <- import("/home/hy395/programs/genomes/blacklists/hg38.blacklist.bed.gz")
}
if (genome == "mm9") {
 blacklist <- import("/home/hy395/programs/genomes/blacklists/mm9.blacklist.bed.gz")
}
if (genome == "mm10") {
 blacklist <- import("/home/hy395/programs/genomes/blacklists/mm10.blacklist.bed.gz")
}
# remove blacklist regions
gr <- gr[countOverlaps (gr, blacklist) == 0]

export.bed(gr, gsub("txt","bed",input))
saveRDS(gr, gsub("txt","rds",input))

