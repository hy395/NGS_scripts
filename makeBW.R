# bam to bigwig
# edited 9.23.2017, the second argument is extsize (2*shift)!!!

library(rtracklayer)
library(GenomicAlignments)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly=TRUE)
print(args)
inputFile <- args[1]
predictd <- as.numeric(args[2])
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


