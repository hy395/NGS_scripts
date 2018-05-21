#merge peaks
collapse.pairwise.celltype.peaks <- function (peaks1, peaks2, overlap.ratio.cutoff=0.75) {
  ## Function to find overlapping ratios
  find.overlap.ratio <- function (gr1, gr2) {
    ## Find Overlaps
    overlaps <- as.matrix (findOverlaps (gr1, gr2))
    if (nrow(overlaps)==0){
      output = NULL
    } else {
      ## Build ranges
      ranges <- cbind (start (gr1)[overlaps[,1]], end (gr1)[overlaps[,1]],
                       start (gr2)[overlaps[,2]], end (gr2)[overlaps[,2]])
      ranges <- t(apply (ranges, 1, sort))
      ## Min widths
      widths <- pmin (width (gr1)[overlaps[,1]], width (gr2)[overlaps[,2]])
      ## Overlap ratios
      overlap.ratio <- (ranges[,3] - ranges[,2])/widths
      #for querys mapping to the same subject, find the rows in overlaps that represent the query with best mapping
      best.gr1 <- tapply (1:nrow (overlaps), overlaps[,1], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
      #for subjects mapping to the same query, find the rows in overlaps that represent the subject with best mapping
      best.gr2 <- tapply (1:nrow (overlaps), overlaps[,2], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
      common <- intersect (best.gr1, best.gr2)      
      output <- list (overlap.ratio = overlap.ratio,
                      ranges = ranges,
                      overlaps = overlaps,
                      best.common = common)
    } 
    return (output)
  }
  ## Reset metadata
  values (peaks1) <- values (peaks2) <- NULL
  ## Overlapping peaks which exceed overlap.ratio.cutoff
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common[or$overlap.ratio[or$best.common] > overlap.ratio.cutoff]
  union <- GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                    IRanges (or$ranges[common,1], or$ranges[common,4])) # modified from Manu's code. He's taking the intersection, I am taking the union.
  ## Overlapping peaks with ratio < overlap.ratio.cutoff
  peaks1 <- peaks1[countOverlaps (peaks1, union) == 0]
  peaks2 <- peaks2[countOverlaps (peaks2, union) == 0]    
  or <- find.overlap.ratio (peaks1, peaks2)
  if (is.null(or)){
    common <- or$best.common
    union <- c(union,GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                              IRanges (or$ranges[common,1], or$ranges[common,4])))    
  }
  ## Non overlapping peaks
  union <- c(union, peaks1[countOverlaps (peaks1, union) == 0])
  union <- c(union, peaks2[countOverlaps (peaks2, union) == 0])
  return (union)
}
mergeCTpeaks <- function(pkLst) {
  combined.peaks <- collapse.pairwise.celltype.peaks (pkLst[[1]], pkLst[[2]])
  if (length(pkLst)>2){
    for (i in 3:length(pkLst)){
      combined.peaks <- collapse.pairwise.celltype.peaks(combined.peaks,pkLst[[i]])
      print(i)
    }
  }
  ## annotate each sample where peak was called.
  for (sample in names(pkLst)) {    
    overlaps <- countOverlaps (combined.peaks, pkLst[[sample]])
    mcols (combined.peaks)[,paste0(sample, ".peak")] <- 0
    mcols (combined.peaks)[,paste0(sample, ".peak")][overlaps > 0] <- 1
  }
  names(combined.peaks) <- 1:length(combined.peaks)
  combined.peaks$pattern <- apply(as.matrix(mcols(combined.peaks)), 1, function(row) {
    argsToPaste0 <- as.list(row)
    do.call(paste0, argsToPaste0)
  })
  return(combined.peaks)
}

############################
# add annotation for atlas #
############################
# atlas needs to be named
library(org.Hs.eg.db)
makeGTF1 <- function(){
  annotationTable<-fread("~/programs/genomes/hg19.refseq.gtf")
  temp<-gsub("gene_id \"","",annotationTable$V9)
  annotationTable$gene_id<-gsub("\".*","",temp)
  annotationTable$symbol<-NA
  temp<-gsub(".*transcript_id \"","",annotationTable$V9)
  annotationTable$transcript_id<-gsub("\";","",temp)
  map <- select(org.Hs.eg.db, keys=annotationTable$gene_id, columns=c("REFSEQ","SYMBOL"), keytype="REFSEQ")
  map <- map[!duplicated(map[,1]),]
  rownames(map) <- map[,1]
  annotationTable$symbol <- map[annotationTable$gene_id,2] 
  a <- annotationTable[,list(seqname=V1,annot=V3,start=V4,end=V5,strand=V7,symbol=symbol,transcript_id=transcript_id)]
  a$width<-abs(a$end-a$start)+1
  a<-a[a$seqname%in%c(paste0("chr", 1:22), "chrX", "chrY"),]
  return(a)
}

makeGTF2<-function(){
  annot1<-readRDS("~/programs/genomes/Homo_sapiens.GRCh37.75.rds")
  output<-data.table(seqname = as.vector(annot1@seqnames), annot = annot1$type, start = start(annot1), end = end(annot1),
                     strand=as.vector(strand(annot1)),symbol=annot1$gene_name,transcript_id=annot1$transcript_id,width=width(annot1))
  output<-output[!is.na(transcript_id),list(seqname,annot,start,end,strand,symbol,transcript_id,width)]
  output<-output[output$seqname%in%c(1:22, "X", "Y"),]
  output$seqname<-paste0("chr",output$seqname)
  mapping<-unique(mcols(annot1)[,c("transcript_id","gene_name","gene_biotype")])
  mapping<-mapping[!is.na(mapping[,1]),]
  rownames(mapping)<-mapping[,1]
  return(list(output,mapping))
}
annotateSites <- function(sitesGR, annotationTab) {  
  # change it to a transcript table, where all rows of the same transcript becomes the same, with the same start and end site
  #transcripts <- annotationTab[, list(start = min(start), end = max(end), symbol = symbol), by = list(transcript_id, seqname, strand)]
  transcripts <- unique(annotationTab[, list(start = min(start), end = max(end), symbol = symbol), by = list(transcript_id, seqname, strand)])
  
  # change it to GR
  transcriptsGR <- GRanges(seqnames = transcripts$seqname,
                           ranges = IRanges(start = transcripts$start,end = transcripts$end),
                           strand = transcripts$strand,
                           transcript_id = transcripts$transcript_id,
                           symbol = transcripts$symbol)
  names(transcriptsGR) <- transcriptsGR$transcript_id
  
  #################
  # Promoter sites.
  #################
  tss <- rbind(transcripts[strand == "+", list(transcript_id, seqname, strand, pos = start, symbol)],
               transcripts[strand == "-", list(transcript_id, seqname, strand, pos = end, symbol)])
  tssGR <- GRanges(seqnames = tss$seqname,
                   ranges = IRanges(start = tss$pos,end = tss$pos),
                   strand = tss$strand,
                   transcript_id = tss$transcript_id,
                   symbol = tss$symbol)
  names(tssGR) <- tssGR$transcript_id
  strand(tssGR) <- "*" # ignore the strand of tssGR
  
  #one index for each row of sitesGR, shows where is its closest tss
  nearestTss <- nearest(x = sitesGR,subject = tssGR)
  distToNearestTss <- distance(x = sitesGR,y = tssGR[nearestTss])
  nearestTssTab <- data.table(peak_id = names(sitesGR),
                              transcript_id = names(tssGR[nearestTss]),
                              dist = distToNearestTss)
  matchTranscript2SiteTab <- nearestTssTab[dist <= 2000, list(peak_id, transcript_id, annot = "promoter")]
  
  ###################
  # Intergenic sites.
  ###################
  # sites that are not in matchTranscript2SiteTab(promoter), and not in transcripts
  strand(transcriptsGR) <- "*" # ignore strand of transcriptsGR from this point.
  overlapWithTranscr <- subsetByOverlaps(sitesGR[!names(sitesGR) %in% matchTranscript2SiteTab$peak_id],transcriptsGR)
  intergPks <- sitesGR[!names(sitesGR) %in% union(matchTranscript2SiteTab$peak_id, names(overlapWithTranscr))]
  nearestTransc <- nearest(x = intergPks, subject = transcriptsGR)
  # add intergenic sites to the table
  matchTranscript2SiteTab <- rbind(matchTranscript2SiteTab, data.table(peak_id = names(intergPks),transcript_id = names(transcriptsGR[nearestTransc]),annot = "intergenic"))
  
  #############
  # Exon sites.
  #############
  exonsGR <- GRanges(seqnames = annotationTab[annot == "exon"]$seqname,
                     ranges = IRanges(start = annotationTab[annot == "exon"]$start,
                                      end = annotationTab[annot == "exon"]$end),
                     strand = "*",
                     transcript_id = annotationTab[annot == "exon"]$transcript_id,
                     symbol = annotationTab[annot == "exon"]$symbol)
  names(exonsGR) <- 1:length(exonsGR)
  ovlWithExon <- findOverlaps(overlapWithTranscr,exonsGR)
  ovlWithExonTab <- data.table(peak_id = names(overlapWithTranscr[queryHits(ovlWithExon)]),
                               exon_id = names(exonsGR[subjectHits(ovlWithExon)]))
  intersectRanges <- pintersect(overlapWithTranscr[ovlWithExonTab$peak_id],
                                exonsGR[ovlWithExonTab$exon_id])
  ovlWithExonTab$intersectWidth <- width(intersectRanges)
  ovlWithExonTab$peak_width <- width(overlapWithTranscr[ovlWithExonTab$peak_id])
  ovlWithExonTab$intersectRatio <- ovlWithExonTab$intersectWidth/ovlWithExonTab$peak_width
  exonPksTab <- ovlWithExonTab[, list(exon_id = exon_id[which.max(intersectRatio)], intersectRatio = max(intersectRatio)), by = peak_id][intersectRatio >= .5]
  exonPksTab$transcript_id <- exonsGR[exonPksTab$exon_id]$transcript_id
  exonPks <- overlapWithTranscr[exonPksTab$peak_id]
  exonPks$transcript_id <- exonPksTab$transcript_id
  matchTranscript2SiteTab <- rbind(matchTranscript2SiteTab,
                                   data.table(peak_id = names(exonPks),transcript_id = exonPks$transcript_id,annot = "exon"))
  
  ###############
  # Intron sites.
  ###############
  intronPks <- sitesGR[setdiff(names(sitesGR), matchTranscript2SiteTab$peak_id)]
  nearestTransc <- nearest(x = intronPks, subject = transcriptsGR)
  matchTranscript2SiteTab <- rbind(matchTranscript2SiteTab,
                                   data.table(peak_id = names(intronPks),transcript_id = names(transcriptsGR[nearestTransc]), annot = "intron"))
  
  return(matchTranscript2SiteTab)
}

dist2TSS <- function(sitesGR, annotationTab) {  
  # change it to a transcript table, where all rows of the same transcript becomes the same, with the same start and end site
  transcripts <- annotationTab[, list(start = min(start), end = max(end), symbol = symbol), by = list(transcript_id, seqname, strand)]
  # change it to GR
  transcriptsGR <- GRanges(seqnames = transcripts$seqname,
                           ranges = IRanges(start = transcripts$start,end = transcripts$end),
                           strand = transcripts$strand,
                           transcript_id = transcripts$transcript_id,
                           symbol = transcripts$symbol)
  names(transcriptsGR) <- transcriptsGR$transcript_id  
  # Promoter sites.
  #################
  # => WARNING: I might be missing some intron guys by not looking specifically
  # at the 5'UTR.
  tss <- rbind(transcripts[strand == "+", list(transcript_id, seqname, strand, pos = start, symbol)],
               transcripts[strand == "-", list(transcript_id, seqname, strand, pos = end, symbol)])
  tssGR <- GRanges(seqnames = tss$seqname,
                   ranges = IRanges(start = tss$pos,end = tss$pos),
                   strand = tss$strand,
                   transcript_id = tss$transcript_id,
                   symbol = tss$symbol)
  names(tssGR) <- tssGR$transcript_id
  
  #one index for each row of sitesGR, shows where is its closest tss
  nearestTss <- nearest(x = sitesGR,subject = tssGR)
  distToNearestTss <- distance(x = sitesGR,y = tssGR[nearestTss])
  nearestTssTab <- data.table(peak_id = names(sitesGR),
                              transcript_id = names(tssGR[nearestTss]),
                              dist = distToNearestTss)
  return(nearestTssTab)
}
add_annot <- function(atlas) {
  # atlas need to be named
  mapping <- makeGTF2()
  match <- annotateSites(atlas, mapping[[1]]) # match is not the same order as atlas
  setkey(match,peak_id);match <- match[names(atlas)] # reorder match
  atlas$annot <- match$annot
  atlas$transcript_id <- match$transcript_id
  atlas$gene_name <- mapping[[2]][match$transcript_id,2]
  atlas$gene_biotype <- mapping[[2]][match$transcript_id,3]
  return(atlas)
}


#count reads
countReads <- function(peaks, readGRLst = NA) {
  if (!is.list(readGRLst)) {
    cat("No reads provided (or not in the form of a list).\n")
    return(0)
  }
  if (is.list(readGRLst)) { 
    for (repName in names(readGRLst)) {
      command <- paste0("peaks$", repName, " <- countOverlaps(peaks, readGRLst[[repName]])")
      eval(parse(text = command))
    } 
  } else {
    cat("Incorrect input format.\n")
    return(0)
  }
  return(peaks)
}

#draw pie chart showing annotation
pieC<-function(v,tf=NULL){
  freq<-table(v)
  png(paste0("figures/",tf,"_pie.png"))
  pie(freq, labels=names(freq), main=paste0("pie chart of annotation of ",tf))
  dev.off()  
}


