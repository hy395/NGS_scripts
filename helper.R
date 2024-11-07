library(parallel)
library(Biostrings)
library(data.table)

##################
# help functions #
##################
get.seqs <- function (org, regions, no.cores = 1) {
  get.seq.chr <- function(chr) {
    seq <- org[[chr]]
    if (class(seq) == "MaskedDNAString") 
      seq <- unmasked(seq)
    chr.regions <- regions[which(as.character(seqnames(regions)) == chr)]
    region.seq <- DNAStringSet(Views(seq, start = start(chr.regions), 
                                     end = end(chr.regions)))
    rev <- as.logical(strand(chr.regions) == "-")
    region.seq[rev] <- reverseComplement(region.seq[rev])
    names(region.seq) <- values(chr.regions)$.ZZ
    gc()
    return(region.seq)
  }
  values(regions)$.ZZ <- sprintf("Unit%d", 1:length(regions))
  seqs <- unique(as.character(seqnames(regions)))
  all.seqs <- mclapply(seqs, get.seq.chr, mc.cores = no.cores)
  all.seqs <- do.call(c, all.seqs)
  inds <- sort(as.numeric(gsub("Unit", "", names(all.seqs))), 
               index.return = TRUE)$ix
  all.seqs <- all.seqs[inds]
  all.seqs <- as.character(all.seqs)
  gc()
  all.seqs <- DNAStringSet(all.seqs)
  gc()
  return(all.seqs)
}

#####################
# find TSS per gene #
#####################
find_tss <- function(gtf) {
    tmp <- gtf[gtf$type=='transcript']

    # if a gene has basic transcripts, consider basic transcripts only
    tmp_dt <- data.table( data.frame(mcols(tmp)))    
    basic_count <- tmp_dt[,sum(tag=='basic'),by=gene_name]
    genes_basic <- basic_count$gene_name[basic_count$V1>1]
    genes_nobasic <- basic_count$gene_name[basic_count$V1==0]
    gr_basic <- tmp[tmp$gene_name %in% genes_basic]
    gr_basic <- gr_basic[gr_basic$tag=='basic']
    gr_nobasic <- tmp[tmp$gene_name %in% genes_nobasic]
    gr <- c(gr_basic, gr_nobasic)

    gr_dt <- data.table(chr=as.character(seqnames(gr)),
                         start=start(gr),
                         end=end(gr),
                         strand=as.character(strand(gr)),
                         gene_name=gr$gene_name)    
    tss_pos <- gr_dt[gr_dt$strand=='+',][, .(pos = min(start)), by = list(gene_name, chr, strand)]
    tss_neg <- gr_dt[gr_dt$strand=='-',][, .(pos = max(end)), by = list(gene_name, chr, strand)]
    tss <- rbind(tss_pos, tss_neg)

    tss_gr <- GRanges(seqnames=tss$chr, ranges=IRanges(tss$pos, tss$pos), strand=tss$strand, gene_name=tss$gene_name)
    return(tss_gr)
}

