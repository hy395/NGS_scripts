# ATAC-seq bam file, shift and smooth, then output .bw

library(rtracklayer)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)
bam <- args[1]
smooth_window <- as.numeric(args[2])

print(paste0("read bam from: ",bam))
print("shift + strand reads by +4, and shift - strand reads by -5")

#read alignment
aln<-readGAlignments(bam)
aln<-as(aln,"GRanges")
aln <- aln[seqnames(aln)%in%paste0("chr",c(1:22,"X","Y"))]
seqlevels(aln) <- paste0("chr",c(1:22,"X","Y"))

#shift alignment by size shift, and take the cut site only
start(aln[strand(aln)=="+"]) <- start(aln[strand(aln)=="+"]) + 4
end  (aln[strand(aln)=="+"]) <- start(aln[strand(aln)=="+"])
end  (aln[strand(aln)=="-"]) <- end(aln[strand(aln)=="-"]) - 5
start(aln[strand(aln)=="-"]) <- end(aln[strand(aln)=="-"])

aln <- resize(aln, fix="center", width = smooth_window)

scalingFactor<-1e6/length(aln)
cov1<-coverage(aln)
cov1<-cov1*scalingFactor
# save the bigwig file
export.bw(cov1, gsub("bam$","bw",bam))
# save the centered and smoothed bam file
#saveRDS(aln, file=gsub("bam$","rds",bam))
