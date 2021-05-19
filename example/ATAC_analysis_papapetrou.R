#############
# functions #
#############
library(rtracklayer)
library(data.table)
library(DESeq2)
source("~/programs/myscripts/ATAC_functions.R")
source("~/programs/myscripts/DESeq2_functions.R")

add_annot <- function(atlas) {
  # use REFSEQ annotation here.
  # only consider annotated genes.
  annot <- makeGTF1()
  annot <- annot[!is.na(symbol),]
  annot <- annot[symbol!=""]
  match <- annotateSites(atlas, annot) # match is not the same order as atlas
  setkey(match,peak_id);match <- match[names(atlas)] # reorder match
  id2symbol <- unique(annot[,list(transcript_id, symbol)])
  setkey(id2symbol, transcript_id)
  match$symbol <- id2symbol[match$transcript_id, symbol]
  atlas$annot <- match$annot
  atlas$transcript_id <- match$transcript_id
  atlas$gene_name <- match$symbol
  return(atlas)
}

##############
# make atlas #
##############
all <- lapply(list.files("../ATAC_idr", pattern="rds", full.names=T), readRDS)
names(all) <- list.files("../ATAC_idr", pattern="rds")
A_combine <- unique(Reduce(union, all[1:3]))
S_combine <- unique(Reduce(union, all[4:6]))
atlas <- collapse.pairwise.celltype.peaks(A_combine, S_combine)
names(atlas) <- paste0("peak", 1:length(atlas))
atlas <- add_annot(atlas)

#########
# DESeq #
#########
all <- lapply(list.files("../ATAC_clean_bam", pattern="rds", full.names=T), readRDS)
names(all) <- gsub(".filter.rmdup.rds","",list.files("../ATAC_clean_bam", pattern="rds"))
all <- countReads(atlas, all)
cnts <- data.frame(mcols(all)[,5:10])
dds <- DESeqDataSetFromMatrix(cnts, 
                              colData=data.frame(cond=factor(c("A","A","A","S","S","S"), levels=c("S","A"))), 
                              design=~cond)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
counts(dds) <- counts(dds) + 1L

# PCA
counts_norm <- counts(dds, normalized=T)
pdf("PCA.pdf", 10, 10)
pca<-prcomp(t(counts_norm))
PCs<-data.frame(pca$x[,1:2])
pVar<-round(pca$sdev^2/sum(pca$sdev^2)*100)
plot(PCs,col=rep(c("red","blue"),each=3),pch=19,main="PCA",xlab=paste0("PC1: ",pVar[1],"% variance"),ylab=paste0("PC2: ",pVar[2],"% variance"))
text(PCs[,1],PCs[,2],labels=rownames(PCs),pos=1)
dev.off()

# differential accessible peaks
dds <- nbinomWaldTest(dds, betaPrior=FALSE)
res <- results(dds,independentFiltering=FALSE)

# scatter
png("scatter.png", 600, 600)
#pdf("scatter.pdf", 7, 7)
x <- rowMeans(counts_norm[,4:6])
y <- rowMeans(counts_norm[,1:3])
label_red <- which(res$log2FoldChange > 1 & res$padj<0.05)
label_blue <- which(res$log2FoldChange < (-1) & res$padj<0.05)
smoothScatter(log2(x), log2(y), colramp = colorRampPalette(c("white", "gray30")),
              xlim=c(-1,12), ylim=c(-1,12), xlab="Adherent", ylab="Suspension", main="Adherent vs. Suspension (ATAC-seq)")
points(log2(x)[label_red], log2(y)[label_red], col="red", pch=20)
points(log2(x)[label_blue], log2(y)[label_blue], col="blue", pch=20)
text(2,15,paste0("up-regulated peaks:",length(label_red),"\ndown-regulated peaks:",length(label_blue)))
text(15,2,"padj < 0.05 \n|log2FC| > 1")
abline(1 ,1, lty=2, col="red")
abline(0,1,col="black")
abline(-1 ,1, lty=2, col="blue")
dev.off()

# table
mcols(atlas) <- data.frame(region=paste0(seqnames(atlas), ":", start(atlas), "-", end(atlas)),
                           data.frame(mcols(atlas)), counts_norm, log2FC=res$log2FoldChange, padj=res$padj)
saveRDS(atlas, file="atlas.rds")
export(atlas, "atlas.bed", format="bed")
write.table(data.frame(mcols(atlas)), file = "atlas.txt")
A_accessible <- data.frame(mcols(atlas[atlas$padj<0.05 & atlas$log2FC>1]))
S_accessible <- data.frame(mcols(atlas[atlas$padj<0.05 & atlas$log2FC<(-1)]))
write.table(S_accessible, file = "S_accessible.txt")
write.table(A_accessible, file = "A_accessible.txt")

# pie chart
colors <- c("#0868ac","#7bccc4","#bae4bc","#f0f9e8")
names(colors) <- c("promoter","intron","exon","intergenic")
toplot1<- table(atlas$annot)
toplot2<- table(atlas[atlas$padj<0.05 & atlas$log2FC>1]$annot)
toplot3<- table(atlas[atlas$padj<0.05 & atlas$log2FC<(-1)]$annot)
pdf("pie.pdf", 20, 6)
par(mfrow=c(1,3))
pie(toplot1, col = colors[names(toplot1)], labels = paste0(names(toplot1),":",toplot1), main = "all peaks")
pie(toplot2, col = colors[names(toplot2)], labels = paste0(names(toplot2),":",toplot2), main = "peaks with higher accessibility in S")
pie(toplot3, col = colors[names(toplot3)], labels = paste0(names(toplot3),":",toplot3), main = "peaks with higher accessibility is A")
dev.off()

# GO analysis
up_genes <- unique(atlas[atlas$padj<0.05 & atlas$log2FC>1]$gene_name)
up_id <- select(org.Hs.eg.db, keys=up_genes, columns="ENTREZID", keytype="SYMBOL")$ENTREZID

down_genes <- unique(atlas[atlas$padj<0.05 & atlas$log2FC<(-1)]$gene_name)
down_id <- select(org.Hs.eg.db, keys=down_genes, columns="ENTREZID", keytype="SYMBOL")$ENTREZID
all_id <- select(org.Hs.eg.db, keys=unique(atlas$gene_name), columns="ENTREZID", keytype="SYMBOL")$ENTREZID

up_res <- mySummary(runGOstats(up_id, all_id), minCount=2)
down_res <- mySummary(runGOstats(down_id, all_id), minCount=2)
write.xlsx(up_res, file="GO_result.xlsx", sheetName="adherent")
write.xlsx(down_res, file="GO_result.xlsx", sheetName="suspension", append=T)


