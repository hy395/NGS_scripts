# find candidate for new atlas
# make plots
library(rtracklayer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(data.table)

# functions
countBin<-function(bam,peakBin){
  vectorized_m<-countOverlaps(peakBin,bam)
  output<-matrix(vectorized_m,ncol=1000,byrow=T)
  return(output)
}
generateMatrix <- function(regions, beds) {  
  peaks_centered_tiles<-tile(regions,1000)
  peaks_centered_tiles_gr<-unlist(peaks_centered_tiles)
  
  # rpm count for each bin
  bin_matrix <- mclapply(beds, function(x) countBin(x,peaks_centered_tiles_gr), mc.cores=40)
  for (i in 1:length(bin_matrix)){
    bin_matrix[[i]]<-bin_matrix[[i]]/sf[i]*1e6
  }
  # put a cap at the intensity
  # quantSig <- quantile(sapply(bin_matrix,function(x) apply(x,1,max)))
  # maxSig <- quantSig["75%"]+1.5*(quantSig["75%"]-quantSig["25%"])
  # bin_matrix <- lapply(bin_matrix, function(x) {
  #   x[x>maxSig] <- maxSig
  #   return(x)
  # })
  bin_matrix <- lapply(bin_matrix, function(x) {
    x[x>1] <- 1
    return(x)
  })
  return(bin_matrix)
}

# reading data
load("../safe_harbor/safe_harbor.rdt")
atlas <- readRDS("analysis/rdata/atac_atlas.rds")
atlas <- atlas[atlas$pattern=="11111"]

names(atlas) <- atlas$name
beds <- lapply(list.files("analysis/rdata", pattern=".rmdup.rds", full.names=T), readRDS)
names(beds) <- gsub("analysis/rdata/","",gsub("_IGO_.*","",list.files("analysis/rdata", pattern=".rmdup.rds", full.names=T)))
beds <- beds[c(1,2,5,3,6,4,7)]

sf <- sapply(beds,length)
hg19 <- GRangesForBSGenome("hg19")
chrSize <- end(ranges(hg19)) ; names(chrSize) <- seqnames(hg19)
chrs <- paste0("chr",c(1:22,"X","Y"))

# atac-seq peaks with safe harbor within 5kb of summit
# 379 left
atlas_10k <- atlas
start(atlas_10k) <- end(atlas_10k) <- start(atlas_10k) + atlas_10k$summit - 1
start(atlas_10k) <- pmax(start(atlas_10k)-5000,1)
end(atlas_10k) <- pmin(end(atlas_10k)+5000 , chrSize[as.character(seqnames(atlas_10k))])
queryHits <- queryHits(findOverlaps(atlas_10k,sf_p_included))
atlas_sf_10k <- atlas_10k[unique(queryHits)]
atlas_sf_original <- atlas[names(atlas_sf_10k)]

# annotation for 379 peaks
#################
# calculate TPM #
#################
# RPM normalized count at peak region
# this is not normalized by peak size
bws <- list.files("analysis/bws", pattern="*.bw", full.names=T)
covs <- mclapply(bws, function(x) import(x,format="bw",as="RleList",which=atlas_sf_original), mc.cores=length(bws))
names(covs) <- gsub("analysis/bws/","",gsub("_IGO_.*","",bws))
maxs <- mclapply(covs, function(x) viewMaxs(as(x[atlas_sf_original], "RleViews")))
maxs <- do.call(cbind, maxs)
mcols(atlas_sf_original) <- maxs
atlas_sf_original <- atlas_sf_original[order(rowMeans(maxs),decreasing=T)]

chain <- import.chain("~/programs/genomes/chains/hg19ToHg38.over.chain")
wins_hg38 <- liftOver(atlas_sf_original, chain)
wins_merged <- list()
for (i in 1:length(wins_hg38)) {
  wins_merged[[i]] <- GRanges(seqnames(wins_hg38[[i]])[1], IRanges(min(start(wins_hg38[[i]])), max(end(wins_hg38[[i]]))), "*")
  print(i)
}
wins_hg38 <- do.call(c, wins_merged)

# generate the annotation table
annotations_379 <- data.frame(peak_id=names(atlas_sf_original), 
           hg19=paste0(seqnames(atlas_sf_original), ":", start(atlas_sf_original), "-", end(atlas_sf_original)),
           hg38=paste0(seqnames(wins_hg38), ":", start(wins_hg38), "-", end(wins_hg38)),
           signal_intensity=rowMeans(data.frame(mcols(atlas_sf_original)))
)
write.csv(annotations_379, file="analysis/output/annotation_379pks.csv")                                                           

################
# plot Heatmap #
################
# tile
# 10k
atlas_sf_10k <- atlas_sf_10k[names(atlas_sf_original)]
bin_matrix_10k <- generateMatrix(atlas_sf_10k, beds)

# indicator where is safe harbor, resolution is 10bp
peaks_centered_tiles<-tile(atlas_sf_10k,1000)
peaks_centered_tiles_gr<-unlist(peaks_centered_tiles)
safe_harbor_m <- countBin(sf_p_included,peaks_centered_tiles_gr)
safe_harbor_m2 <- countBin(sf_p_excluded,peaks_centered_tiles_gr)
p_include_safe <- round(sum(safe_harbor_m==1)/length(safe_harbor_m),2)
p_exclude_safe <- round(sum(safe_harbor_m2==1)/length(safe_harbor_m2),2)

annotation <- cbind(vista=as.numeric(atlas_sf_10k$vista), H3K4me1=as.numeric(atlas_sf_10k$H3K4me1))

# rpm count for each bin
bin_matrix <- mclapply(beds, function(x) countBin(x,peaks_centered_tiles_gr), mc.cores=40)

# plot 1
colors <- colorRamp2(c(0, 0.5, 1), c("blue","white","red"))
ht_list = Heatmap(bin_matrix_10k[[1]], name = names(beds)[1], column_title = names(beds)[1], cluster_rows = FALSE, cluster_columns = FALSE, col = colors)
for (i in 2:length(bin_matrix_10k)){
  temp = Heatmap(bin_matrix_10k[[i]], name = names(beds)[i], column_title = names(beds)[i], cluster_rows = FALSE, cluster_columns = FALSE, col = colors)
  ht_list = ht_list + temp
}
ht_list = ht_list + Heatmap(safe_harbor_m, name = "safe_regions", column_title = paste0("safe_regions:", p_include_safe), 
                            cluster_rows = FALSE,cluster_columns = FALSE, col = colorRamp2(c(0,1), c("blue","green")))
#ht_list = ht_list + Heatmap(safe_harbor_m2, name = "safe_regions_no_p", column_title = paste0("safe_regions_no_p:", p_exclude_safe), 
#                            cluster_rows = FALSE, cluster_columns = FALSE, col = colorRamp2(c(0,1), c("blue","green")))
#ht_list = ht_list + Heatmap(annotation, name = "annotation", column_title = "annotation", 
#                            cluster_rows = FALSE, cluster_columns = FALSE, show_column_names=T, width = unit(2, "cm"), col = colorRamp2(c(0,1), c("white","green")))
png("output/379_atac_atlas.png",1500, 800)
draw(ht_list)
dev.off()

# pks need to > 400bp, so we can acommondate at least 3 guides.
# 293 left
atlas_sf_original <- atlas_sf_original[width(atlas_sf_original) > 400]
saveRDS(atlas_sf_original, file="rdata/atac_safe.rds")

################
# plot Heatmap #
################
# tile
# 10k
atlas_sf_10k <- atlas_sf_10k[names(atlas_sf_original)]
bin_matrix_10k <- generateMatrix(atlas_sf_10k, beds)

# indicator where is safe harbor, resolution is 10bp
peaks_centered_tiles<-tile(atlas_sf_10k,1000)
peaks_centered_tiles_gr<-unlist(peaks_centered_tiles)
safe_harbor_m <- countBin(sf_p_included,peaks_centered_tiles_gr)
safe_harbor_m2 <- countBin(sf_p_excluded,peaks_centered_tiles_gr)
p_include_safe <- round(sum(safe_harbor_m==1)/length(safe_harbor_m),2)
p_exclude_safe <- round(sum(safe_harbor_m2==1)/length(safe_harbor_m2),2)

annotation <- cbind(vista=as.numeric(atlas_sf_10k$vista), H3K4me1=as.numeric(atlas_sf_10k$H3K4me1))


# rpm count for each bin
bin_matrix <- mclapply(beds, function(x) countBin(x,peaks_centered_tiles_gr), mc.cores=40)

# plot 1
colors <- colorRamp2(c(0, 0.5, 1), c("blue","white","red"))
ht_list = Heatmap(bin_matrix_10k[[1]], name = names(beds)[1], column_title = names(beds)[1], cluster_rows = FALSE, cluster_columns = FALSE, col = colors)
for (i in 2:length(bin_matrix_10k)){
  temp = Heatmap(bin_matrix_10k[[i]], name = names(beds)[i], column_title = names(beds)[i], cluster_rows = FALSE, cluster_columns = FALSE, col = colors)
  ht_list = ht_list + temp
}
ht_list = ht_list + Heatmap(safe_harbor_m, name = "safe_regions", column_title = paste0("safe_regions:", p_include_safe), 
                            cluster_rows = FALSE,cluster_columns = FALSE, col = colorRamp2(c(0,1), c("blue","green")))
ht_list = ht_list + Heatmap(safe_harbor_m2, name = "safe_regions_no_p", column_title = paste0("safe_regions_no_p:", p_exclude_safe), 
                            cluster_rows = FALSE, cluster_columns = FALSE, col = colorRamp2(c(0,1), c("blue","green")))
ht_list = ht_list + Heatmap(annotation, name = "annotation", column_title = "annotation", 
                            cluster_rows = FALSE, cluster_columns = FALSE, show_column_names=T, width = unit(2, "cm"), col = colorRamp2(c(0,1), c("white","green")))
png("output/293_atac_atlas.png",2000,800)
draw(ht_list)
dev.off()

########
# additional repeat annotation
tmp <- read.table("external_source/hg19_repeats.txt", comment.char="#", sep="\t", colClasses="character")
repeats <- tmp[grep("LINE|SINE|LTR",tmp$V12),]

repeats <- GRanges(repeats$V6, IRanges(start=pmin(as.numeric(repeats$V7), as.numeric(repeats$V8)),
                            end=pmax(as.numeric(repeats$V7), as.numeric(repeats$V8))),"*")
atlas_10k_tile <- unlist(tile(atlas_sf_10k,1000))
repeat_annot <- rep(0, length(atlas_10k_tile))
repeat_id <- unique(queryHits(findOverlaps(atlas_10k_tile, repeats)))
repeat_annot[repeat_id] <- 1
repeat_m <- matrix(repeat_annot,ncol=1000,byrow=T)
ht_list = ht_list + Heatmap(repeat_m, name = "repeat", column_title = "repeat", 
                            cluster_rows = FALSE, cluster_columns = FALSE, show_column_names=F, width = unit(2, "cm"), col = colorRamp2(c(0,1), c("white","green")))
png("output/293_atac_atlas_2.png",2000,800)
draw(ht_list)
dev.off()
##########


# 
# ###############################
# # make a list of safe harbors #
# ###############################
# # only consider safe harbors that are <5kb from a atac-seq peaks
# # bin those regions into 500 bp bins
# # columns are: 
# # 1. safe harbor genomic location
# # 2. safe harbor width
# # 3. safe harbor average accessibility (ATAC-seq RPM)
# # 4. nearest ATAC-seq id
# # 5. distance to nearest ATAC-seq center
# # 6. width of ATAC-seq peak
# # 7. maximum signal at nearest ATAC-seq
# # 8. distance to the closest transcript
# # 9. name of the closest transcript
# # 10. type of the closest transcript
# readTranscripts <- function(){
#   # read gencode
#   annotationTable <- import("../safe_harbor/raw_data/gencode.v25lift37.annotation.gtf")
#   annotationTable <- annotationTable[annotationTable$type=="gene"]
#   gene <- annotationTable$gene_name;type <- annotationTable$gene_type
#   gencode <- granges(annotationTable); gencode$gene <- gene; gencode$type <- type
#   gencode <- keepSeqlevels(gencode, chrs)
#   # read refseq
#   refseq <- import("../safe_harbor/raw_data/hg19_refseq.gtf.gz",format="gtf")
#   refseq_dt <- data.table(seqnames=as.character(seqnames(refseq)), start=start(refseq), end=end(refseq), strand=as.character(strand(refseq)), gene_id=refseq$gene_id, transcript_id=refseq$transcript_id)
#   refseq_dt_by_tran <- unique(refseq_dt[, list(start = min(start), end = max(end), gene_id=gene_id), by = list(transcript_id, seqnames, strand)])
#   refseq <- GRanges(seqnames = refseq_dt_by_tran$seqnames,
#                     ranges = IRanges(start = refseq_dt_by_tran$start,end = refseq_dt_by_tran$end),
#                     strand = refseq_dt_by_tran$strand,
#                     gene = refseq_dt_by_tran$gene_id)
#   refseq <- unique(sort(refseq))
#   refseq <- refseq[grep("NM_", refseq$gene)]
#   refseq$type <- gsub("_.*","",refseq$gene)
#   refseq <- keepSeqlevels(refseq, chrs)
#   all <- c(gencode,refseq)
#   return(all)
# }
# all <- readTranscripts()
# bws <- list.files("bws", pattern="*.bw", full.names=T)
# covs <- mclapply(bws, function(x) import(x,format="bw",as="RleList"), mc.cores=length(bws))
# names(covs) <- gsub("bws/","",gsub("_IGO_.*","",bws))
# 
# tmp <- sf_p_excluded
# # 1. safe harbor
# sf1 <- intersect(tmp, atlas_10k)
# sf1 <- sort(c(sf1[width(sf1)<600], unlist(tile(sf1[width(sf1)>=600],width=500))))
# pos <- paste0(seqnames(sf1),":", start(sf1),"-", end(sf1)) # column 1
# width_sf <- width(sf1) # column 2
# # safe harbor average rpm
# means <- mclapply(covs, function(x) viewMeans(as(x[sf1], "RleViews")))
# sf_mean_rpm <- rowMeans(do.call(cbind,means)) # column 3
# # distance to nearest atac-seq center
# centers <- resize(atlas,fix="center",width=1)
# dist_to_atac <- distanceToNearest(sf1, centers, ignore.strand=TRUE)
# atac_id <- names(atlas[dist_to_atac@subjectHits]) # column 4
# distance_to_neareast_atac_center <- mcols(dist_to_atac)$distance # column 5
# # annotation for the nearest atac peak (ovlp with vista and ovlp with H3K4me1)
# vista <- atlas[atac_id]$vista # column 6
# h3k4me1 <- atlas[atac_id]$H3K4me1 # column 7
# # nearest atac-seq width
# atac_width <- width(atlas[dist_to_atac@subjectHits]) # column 8
# # nearest atac-seq max coverage (across 6 cell types)
# # two atlas appear out of the coverage, first remove those two and finally add their coverage back as 0
# bad <- unique(unlist(lapply(covs, function(x) {
#   a <- sapply(x,length)
#   b <- findOverlaps(atlas,GRanges(seqnames=names(a), ranges=IRanges(1,a), strand="*"), type="within")
#   return(which(!c(1:length(atlas))%in%b@queryHits))
# })))
# a <- atlas
# if(length(bad)>0) {
#   a <- atlas[-bad]  
# }
# maxs <- mclapply(covs, function(x) viewMaxs(as(x[a], "RleViews")))
# all_max <- rowMeans(do.call(cbind, maxs))
# # all_max <- rep(0,length(atlas))
# # all_max[as.numeric(names(a))] <- maxs
# names(all_max) <- names(atlas)
# max_sig_of_nearest_atac <- all_max[dist_to_atac@subjectHits] # column 9
# # distance to closest transcript
# dist_to_tran <- distanceToNearest(sf1, all, ignore.strand=TRUE)
# dist_to_nearest_tran <- mcols(dist_to_tran)$distance # column 10
# # transcript annotation
# trans <- data.frame(mcols(all[dist_to_tran@subjectHits])) # column 11 and 12
# 
# ########################
# # make the final table #
# ########################
# final <- data.frame(sf_pos=pos, sf_size=width_sf, sf_mean_rpm, 
#                     atac_id=atac_id, vista=vista, h3k4me1=h3k4me1, atac_dist=distance_to_neareast_atac_center, atac_width, atac_sig=max_sig_of_nearest_atac, 
#                     gene_dist=dist_to_nearest_tran, gene_name=trans$gene, gene_type=trans$type)
# final <- final[order(final$atac_sig,decreasing = T),]
# saveRDS(final, "output/final_candidates_exclude_p.rds")
# write.csv(final, file="output/final_candidates_exclude_p.csv")
