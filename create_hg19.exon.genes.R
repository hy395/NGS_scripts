library(rtracklayer)

#hg19
chrs<-c(1:22,"X","Y")
#annotation GRList
annotations<-import('Homo_sapiens.GRCh37.75.gtf.gz')
annotations<-keepSeqlevels(annotations,chrs)
genes <- unique(annotations$gene_name)
exons <- annotations[annotations$type == 'exon']
exons.genes <- disjoin(split(exons, factor(exons$gene_name,levels=genes)))
tRNA_rRNA <- import("hg19_trna_rrna.gtf", format="gtf",asRangedData = FALSE)
seqlevels(tRNA_rRNA)<-gsub("chr","",seqlevels(tRNA_rRNA))
tRNA_rRNA<-keepSeqlevels(tRNA_rRNA,chrs)
rrnas <-unique(tRNA_rRNA$gene_id)
tRNA_rRNA <- disjoin(split(tRNA_rRNA, factor(tRNA_rRNA$gene_id,levels=rrnas)))
tRNA_rRNA <- endoapply(tRNA_rRNA, function (x) {mcols(x) <- NULL; x})
overlap<-findOverlaps(exons.genes,tRNA_rRNA)
rrnas<-unique(overlap@queryHits)
exons.genes<-exons.genes[-rrnas]
saveRDS(exons.genes, file="hg19.exons.genes.rds")

