# analysis
library(DESeq2)
library(GenomicAlignments)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
source("~/programs/myscripts/DESeq2_functions.R")

#########
# input #
#########
counts <- readRDS("counts.rds")
counts <- updateObject(counts)
colnames(counts) <- gsub("_R1.filter.bam","",colnames(counts))
cnts <- assay(counts)
print(colSums(cnts))

# only look at coding genes
human <- useMart(host="www.ensembl.org",biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl")
map <- getBM(attributes=c("hgnc_symbol","gene_biotype"), filters ="hgnc_symbol",values=rownames(counts), mart=human)
map <- map[map$gene_biotype=="protein_coding",]
counts <- counts[rownames(counts)%in%map[,1], ]

#look at genes with highest counts to see if we have any tRNA or rRNA left
a<-assay(counts)
a[order(apply(a,1,max),decreasing=T)[1:20],]

#######################
# create DESeq object #
#######################
colData(counts)$cond <- factor(c(rep("A",3),rep("S",3)), levels=c("S","A"))
dds<-DESeqDataSet(counts,design=~cond)
assays(dds)$fpkm <- fpkm(dds)

# rpkm filter
fpkmFilterRange(dds,minSample=3,fpkmRange=c(0,0.1,0.5,1))
dds_filter <- fpkmFilter(dds,minSample=3,minFPKM=0.5)

# count filter
dds_filter <- estimateSizeFactors(dds_filter)
countFilterRange(dds_filter,minSample=3,countRange=c(0,1,10,100))
dds_filter <- countFilter(dds_filter,minSample=3,minCount=10)
dds_filter <- estimateSizeFactors(dds_filter)
dds_filter <- estimateDispersions(dds_filter)
counts(dds_filter) <- counts(dds_filter) + 1L

############
# draw PCA #
############
counts_norm <- counts(dds_filter, normalized=T)
pdf("PCA.pdf", 10, 10)
pca<-prcomp(t(counts_norm))
PCs<-data.frame(pca$x[,1:2])
pVar<-round(pca$sdev^2/sum(pca$sdev^2)*100)
plot(PCs,col=rep(c("red","blue"),each=3),pch=19,main="PCA",xlab=paste0("PC1: ",pVar[1],"% variance"),ylab=paste0("PC2: ",pVar[2],"% variance"))
text(PCs[,1],PCs[,2],labels=rownames(PCs),pos=1)
dev.off()

##################
# DESeq analysis #
##################
dds_filter <- nbinomWaldTest(dds_filter, betaPrior=FALSE)
res <- results(dds_filter,independentFiltering=FALSE)
counts_norm <- counts(dds_filter, normalized=T)
saveRDS(dds_filter, file="dds.rds")

################
# plot scatter #
################
pdf("scatter.pdf", 8, 8)
x <- rowMeans(counts_norm[,4:6])
y <- rowMeans(counts_norm[,1:3])
label_red <- which(res$log2FoldChange > 0 & res$padj<0.05)
label_blue <- which(res$log2FoldChange < 0 & res$padj<0.05)
smoothScatter(log2(x), log2(y), colramp = colorRampPalette(c("white", "gray30")),
              xlim=c(-1,17), ylim=c(-1,17), xlab="Suspension", ylab="Adherent", main="Adherent vs. Suspension (RNA-seq)")
points(log2(x)[label_red], log2(y)[label_red], col="red", pch=20)
points(log2(x)[label_blue], log2(y)[label_blue], col="blue", pch=20)
text(2,15,paste0("up-regulated genes:",length(label_red),"\ndown-regulated genes:",length(label_blue)))
text(15,2,"padj < 0.05")
abline(1 ,1, lty=2, col="red")
abline(0,1,col="black")
abline(-1 ,1, lty=2, col="blue")
dev.off()

##########
# GOstat #
##########
sus <- rownames(res)[res$log2FoldChange > 0 & res$padj<0.05]
adh <- rownames(res)[res$log2FoldChange < 0 & res$padj<0.05]
gene_id <- select(org.Hs.eg.db, keys = adh, columns = "ENTREZID", keytype = "SYMBOL")
all_id <- select(org.Hs.eg.db, keys = rownames(res), columns = "ENTREZID", keytype = "SYMBOL")
gostat_res <- mySummary(runGOstats(unique(gene_id[,2]), unique(all_id[,2])))

#################
# save to table #
#################
# log2FC > 1, padj < 0.05
towrite <- data.frame(counts_norm, data.frame(res)[,c(2,6)])
up <- towrite[label_red,]; up <- up[order(up$padj), ]
down <- towrite[label_blue,]; down <- down[order(down$padj), ]
write.table(up, file="A.txt")
write.table(down, file="S.txt")
write.table(towrite, file="all.txt")

##################
# Heatmap of DEG #
##################
# log2FC > 2, padj < 0.05
padj_thres <- 0.05
log2fc_thres <- 0
up <- rownames(res[which((res$padj<padj_thres) & (res$log2FoldChange>log2fc_thres)),])
down <- rownames(res[which((res$padj<padj_thres) & (res$log2FoldChange<(-log2fc_thres))),])
toplot <- counts_norm[c(up,down),c(4:6,1:3)]
labels <- as.character(colData(dds)$cond)[c(4:6,1:3)]
ha=HeatmapAnnotation(df = data.frame(type=as.character(labels)), col = list(type = c("S"="#cbd5e8","A" = "#fdcdac")))
#pdf("Heatmap1_3.pdf", 14, 10)
#ht1 <- Heatmap(t(scale(t(log(toplot[1:34,])))), show_row_names=T, top_annotation=ha, show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 10), cluster_rows=T, column_names_max_height=unit(6, "cm"))
#ht2 <- Heatmap(t(scale(t(log(toplot[35:68,])))), show_row_names=T, top_annotation=ha, show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 10), cluster_rows=T, column_names_max_height=unit(6, "cm"))
#ht3 <- Heatmap(t(scale(t(log(toplot[69:102,])))), show_row_names=T, top_annotation=ha, show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 10), cluster_rows=T, column_names_max_height=unit(6, "cm"))
#ht1+ht2+ht3
#dev.off()
#pdf("Heatmap4.pdf", 6, 10)
#ht4 <- Heatmap(t(scale(t(log(toplot[103:137,])))), show_row_names=T, top_annotation=ha, show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 10), cluster_rows=T ,column_names_max_height=unit(6, "cm"))
#ht4
#dev.off()
pdf("Heatmap_all.pdf", 8, 8)
Heatmap(t(scale(t(log(toplot)))), show_row_names=T, top_annotation=ha, show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 0.3), cluster_rows=T ,column_names_max_height=unit(6, "cm"))
dev.off()


##########
# GOstat #
##########


########
# GSEA #
########
genGseaFiles(dds_filter, colData(dds_filter)$cond, "GSEA/")
runGSEA("GSEA/labels.cls", "GSEA/expr.txt", "~/programs/genomes/msigdb/c2.all.v6.1.symbols.gmt", "#A_versus_S", 1024, "gene_set", "c2", "GSEA/c2")
runGSEA("GSEA/labels.cls", "GSEA/expr.txt", "~/programs/genomes/msigdb/c5.all.v6.1.symbols.gmt", "#A_versus_S", 1024, "gene_set", "c5", "GSEA/c5")
runGSEA("GSEA/labels.cls", "GSEA/expr.txt", "~/programs/genomes/msigdb/c7.all.v6.1.symbols.gmt", "#A_versus_S", 1024, "gene_set", "c7", "GSEA/c7")

