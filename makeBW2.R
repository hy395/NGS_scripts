# this script takes a bam file, and pred_d file.
# and outputs the bigwig file.

library(rtracklayer)
library(GenomicAlignments)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly=TRUE)
print(args)

inputFile <- args[1]
inputPred <- args[2]

# extract predictd
pred <- read.table(inputPred, sep="\t", comment.char="", stringsAsFactors = F)
tmp <- pred$V1[grep("predicted fragment length",pred$V1)]
d <- gsub(" bps ","",gsub(".*predicted fragment length is ","",tmp))
predictd <- as.numeric(d)
if(is.na(predictd)) {
  stop("fail to extract shift size")
} else {
  print(paste0("extracted d is: ", predictd))
}

shift <- floor(predictd/2)
output <- gsub("bam$","bw",inputFile)

readBam<-function(bam,shift,bigWig=NULL){
  print(paste0("read bam from: ",bam))
  print(paste0("shift length is: ",shift))
  print(paste0("save output bigwig to: ",bigWig))
  
  #specify an output file name if want to output a BigWig coverage file.
  #read alignment
  aln<-readGAlignments(bam)
  aln<-as(aln,"GRanges")
  
  #shift alignment by size shift
  aln[strand(aln)=="+"]<-shift(aln[strand(aln)=="+"],shift)
  aln[strand(aln)=="-"]<-shift(aln[strand(aln)=="-"],-shift)
  
  #scale to get read per million and output to BigWig
  if(!is.null(bigWig)){
    scalingFactor<-1e6/length(aln)
    #calculate coverage
    cov1<-coverage(aln)
    cov1<-cov1*scalingFactor
    export.bw(cov1,bigWig)    
  }  
  return(aln)
}

temp <- readBam(inputFile,shift,output)


