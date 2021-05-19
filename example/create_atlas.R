#########################
# merge to create atlas #
#########################

library(rtracklayer)

filter.idr <- function(input, p.value) {
  a <- read.table(input, sep="\t", header=F, stringsAsFactors=F, colClasses="character")
  gr <- GRanges(seqnames=a$V1, ranges=IRanges(as.numeric(a$V2), as.numeric(a$V3)), strand="*", 
                score=as.numeric(a$V5), signalValue=as.numeric(a$V7), summit=as.numeric(a$V10), localIDR=as.numeric(a$V11), globalIDR=as.numeric(a$V12))
  thres <- floor(-125*log2(p.value))
  gr <- gr[seqnames(gr)%in%paste0("chr",c(1:22,"X","Y"))]
  gr <- gr[gr$score>=thres]
  
  return(gr)  
}

# create the atlas, peaks with that is reproducible in at least one idr
files <- list.files("idr", pattern="txt", full.names=T)
pks <- lapply(files, function (x) {
  return(filter.idr(x, 0.05))
})
atlas <- import("idr/all_pool_peaks.narrowPeak",format="BED",extraCols = c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", summit = "integer"))
ovlps <- lapply(pks, function(x) {
  a <- rep(0, length(atlas))
  a[unique(findOverlaps(atlas, x)@queryHits)] <- 1
  return(a)
})
cell.types <- do.call(cbind, ovlps)
colnames(cell.types) <- c("CD3.1","CD3.2","CD3.3","CD4","CD8")
mcols(atlas) <- data.frame(mcols(atlas)[,c(1,6)], cell.types)
atlas_filter <- atlas[rowSums(data.frame(mcols(atlas)[,3:7]))!=0,]
atlas_filter$pattern <- do.call(paste0,data.frame(mcols(atlas_filter)[,3:7]))
mcols(atlas_filter) <- mcols(atlas_filter)[,-c(3:7)]

# remove blacklist regions
blacklist <- import("~/analysis/genomes/wgEncodeDacMapabilityConsensusExcludable.bed")
atlas_filter <- atlas_filter[countOverlaps (atlas_filter, blacklist) == 0]

# annotate ovlp with vistaEnhancer
load("external_source/external.rdt")
atlas_filter$vista <- "0"
atlas_filter$vista [countOverlaps (atlas_filter, enhancer) != 0] <- "neg"
atlas_filter$vista [countOverlaps (atlas_filter, enhancer[enhancer$annotation=="positive"]) != 0] <- "pos"

# annotate ovlp with ENCODE CD34 H3K4me1
atlas_filter$H3K4me1 <- countOverlaps (atlas_filter, CD3_H3K4me1) != 0

# output to bed file, output only peaks that are called in ALL IDRs
output <- atlas_filter[atlas_filter$pattern=="11111"]

export.bed(atlas_filter, "output/pks_atlas.bed")
export.bed(output, "output/pks_all_rep.bed")

saveRDS(atlas_filter, "rdata/atac_atlas.rds")

