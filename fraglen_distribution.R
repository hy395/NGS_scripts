################################
# compute fragment length dist #
################################
# run on _SE.bed files
library(data.table)

options(echo=TRUE)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>1) stop("only takes single input file")
outfile <- gsub("\\.bed" ,"_fraglen.pdf", args)

# read data
system(paste0("cat ",args," | cut -f2-3 | sed 'N;s/\\n/\\t/' > tmp.txt"))
tmp <- fread("tmp.txt")
system("rm tmp.txt")

# compute fragment length
tmp[,min_val := pmin(V1,V2,V3,V4)]
tmp[,max_val := pmax(V1,V2,V3,V4)]

fraglen <- tmp$max_val - tmp$min_val
fraglen <- fraglen[fraglen<1000]

pdf(outfile)
plot(density(fraglen), main="fragment length distribution", xlab="fragment length")
dev.off()

