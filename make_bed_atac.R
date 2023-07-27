# ATAC-seq bam file, shift and smooth, then output .bw

library(rtracklayer)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)
bam <- args[1]

print(paste0("read bam from: ",bam))
print("shift + strand reads by +4, and shift - strand reads by -5")

#read alignment
aln<-readGAlignments(bam)
print("finish reading bam file.")

aln<-as(aln,"GRanges")
aln <- aln[seqnames(aln)%in%paste0("chr",c(1:22,"X","Y"))]

#shift alignment by size shift, and take the cut site only
start(aln[strand(aln)=="+"]) <- start(aln[strand(aln)=="+"]) + 4
end  (aln[strand(aln)=="+"]) <- start(aln[strand(aln)=="+"])
end  (aln[strand(aln)=="-"]) <- end(aln[strand(aln)=="-"]) - 5
start(aln[strand(aln)=="-"]) <- end(aln[strand(aln)=="-"])

export.bed(aln, gsub("bam$","bed", bam))

#print("shifted bam file.")

#scalingFactor<-length(aln)/1e6
#cov1<-coverage(aln)
#cov1<-cov1/scalingFactor
#print("finish creating coverage file.")

# save the bigwig file
#export.bedgraph(cov1, gsub("bam$","bedgraph",bam))
