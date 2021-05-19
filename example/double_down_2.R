#################
# 42 candidates #
#################
library(data.table); library(rtracklayer); library(DESeq2); library(dplyr);library(org.Hs.eg.db);library(ggplot2);library(xlsx);library(gtools)
source("~/programs/myscripts/ATAC_functions.R")
source("~/programs/myscripts/helper.R")

atlas <- readRDS("../ATAC_analysis/atlas.rds")
# read a annotation table
annot <- makeGTF1()
annot <- annot[!is.na(symbol)]
annot <- annot[symbol!=""]
# remove those with no gene symbol
genes <- unique( annot[, list(start = min(start), end = max(end), symbol = symbol), by = list(symbol, seqname, strand)])
genesGR <- GRanges(seqnames = genes$seqname,
                   ranges = IRanges(start = genes$start,end = genes$end),
                   strand = genes$strand,
                   symbol = genes$symbol)
names(genesGR) <- genesGR$symbol
nearestGene <- nearest(x = atlas,subject = genesGR, ignore.strand=T)
atlas$gene_name <- genesGR[nearestGene]$symbol
atlas$distToNearest <- distance(x = atlas,y = genesGR[nearestGene])

# make RNA dataframe
rna <- readRDS("../RNA_analysis/dds.rds")
res <- results(rna,independentFiltering=FALSE)
counts_norm <- counts(rna, normalized=T)
rna <- data.frame(counts_norm, padj=res$padj, log2FC=res$log2FoldChange)
rna <- rna[order(rna$log2FC, decreasing=T), ]
rna$gene <- rownames(rna)

# make ATAC dataframe
peaks <- data.frame(mcols(atlas))
peaks$name <- names(atlas)
peaks$gene <- peaks$gene_name
peaks$gene.dist <- atlas$distToNearest
peaks <- peaks[peaks$distToNearest < 50000, ]
#peaks <- peaks[peaks$annot == "promoter", ]
median_ATAC_log2FC <- tapply(peaks$log2FC, peaks$gene, median)

###############
# scatterplot #
###############
ovlp <- intersect(rownames(rna), names(median_ATAC_log2FC))
x <- rna[ovlp, "log2FC"]
y <- median_ATAC_log2FC[ovlp]
pdf("double_down_2/rna_atac_scatter_promoter.pdf")
plot(x, y, pch=20)
text(2,2,paste0("coorelation is: ", round(cor(x, y),3)))
dev.off()

#################################
# KS test for identifying genes #
#################################
# how many genes have only 1 peak? 6816
pks_by_gene <- split(peaks, f=peaks$gene)
res <- lapply(pks_by_gene, function(x){
  a <- ks.test(x$log2FC, peaks$log2FC)
  return(c(a$statistic, a$p.value))
})
saveRDS(res, file="double_down_2/gene_level_accessibility_KS.rds")

#################################
# focus on genes up in adherent #
#################################
# that is, the double down genes

rna_up <- rownames(rna)[rna$padj < 0.05 & rna$log2FC > 0]
rna_down <- rownames(rna)[rna$padj < 0.05 & rna$log2FC < 0]

peak_up1 <- peaks[peaks$padj<0.05 & peaks$log2FC>1, ]$gene_name
peak_up2 <- intersect(names(median_ATAC_log2FC[median_ATAC_log2FC>0]), rownames(res)[res[,2]<0.05])
peak_down1 <- peaks[peaks$padj<0.05 & peaks$log2FC<(-1), ]$gene_name
peak_down2 <- intersect(names(median_ATAC_log2FC[median_ATAC_log2FC<0]), rownames(res)[res[,2]<0.05])

# double_up
peak_up <- union(peak_up1, peak_up2)
double_up <- intersect(rna_up, peak_up)
tmp <- read.xlsx("double_down_2/double_down_genes_old.xlsx", sheetIndex=1)
toremove <- double_up[!double_up %in% as.character(tmp[-1,2])]
peak_up <- peak_up[!peak_up %in% toremove]
double_up <- intersect(rna_up, peak_up)

peak_down <- union(peak_down1, peak_down2)
double_down <- intersect(rna_down, peak_down)

draw_diamond <- function(rna.data, atac.data, genes=NULL, top.n=NULL, max.gene.dist=50000, diamond.dist=0.1) {
  # inspired by Yuri's code for drawing
  # @rna.data requires column: gene, log2FC, padj
  # @atac.data requires column: gene, gene.dist, log2FC, padj
  # @either genes or top.n need to be specified
  
  # data preparation
  rna.data <- rna.data[rna.data$gene %in% atac.data$gene,] # only consider genes with atac peak
  rownames(rna.data) <- rna.data$gene
  atac.data <- atac.data[atac.data$gene.dist <= max.gene.dist, ]
  
  # plot1 is RNA-seq for selected genes
  if (!is.null(genes)) {
    toplot1 <- rna.data[genes, c("gene","log2FC","padj")]
    toplot1 <- toplot1[order(toplot1$log2FC), ]
    toplot1$rank <- 1:nrow(toplot1)
  }
  if (!is.null(top.n)) {
    toplot1 <- rna.data[, c("gene","log2FC","padj")]
    toplot1 <- toplot1[order(toplot1$log2FC), ]
    toplot1 <- toplot1[c(1:top.n, (nrow(toplot1)-top.n+1):nrow(toplot1)), ]
    toplot1$rank <- 1:nrow(toplot1)
  }
  
  # plot2 is ATAC
  toplot2 <- atac.data[atac.data$gene %in% toplot1$gene, c("gene", "log2FC", "padj")]
  toplot2$x <- toplot1[toplot2$gene,"rank"]
  toplot2$y <- 0
  toplot2 <- split(toplot2, toplot2$gene)
  toplot2 <- lapply(toplot2, function(tmp) {
    tmp <- tmp[order(tmp$log2FC), ]
    y <- rep(toplot1[unique(tmp$gene), "log2FC"], nrow(tmp))
    y <- y + seq(diamond.dist, diamond.dist * nrow(tmp), by = diamond.dist)
    tmp$y <- y
    return(tmp)
  })
  toplot2 <- do.call(rbind, toplot2)
  
  # plot all
  y_max <- max(abs(toplot1$log2FC))+1
  p <- ggplot(toplot1,  aes(rank, log2FC, label = gene)) + 
    geom_text(aes(angle=90), hjust = 1) +
    geom_point(data=toplot2, aes(x, y, color=log2FC)) +
    geom_point(data=toplot2, aes(x, y), shape=1, stroke=0.2) +
    labs(colour = "ATAC_log2FC") +
    scale_colour_gradient2(low="royalblue3", mid="white", high="firebrick3", midpoint = 0) +
    theme_classic() +  ylim(-y_max,y_max) + xlab("") + ylab("RNA_log2FC") +
    geom_hline(yintercept=0, linetype = 2)
  return(p)
  
}

peaks$log2FC [peaks$log2FC>1] <- 1
peaks$log2FC [peaks$log2FC<(-1)] <- (-1)

top <- draw_diamond (rna, peaks, genes=c(double_down, double_up))
pdf("double_down_2/diamond_double_v4.pdf", 18, 8)
print(top)
dev.off()

################
# create table #
################
peaks.dt <- data.table(peaks)
peaks_by_gene <- data.frame(peaks.dt[,.SD[which(log2FC==max(log2FC))],by=gene_name,.SDcols=c("region","annot","distToNearest","padj","log2FC")])
rownames(peaks_by_gene) <- peaks_by_gene$gene_name
RNA_info <- rna[double_up,c("log2FC","padj")]
ATAC_info <- peaks_by_gene[double_up, c("region","annot","distToNearest","log2FC","padj")]
ATAC_KS <- res[double_up, ]
double_up_info <- cbind(RNA_info, ATAC_info, ATAC_KS)
full_name <- select(org.Hs.eg.db, double_up, c("SYMBOL","GENENAME"), keytype="SYMBOL")
double_up_info <- cbind(full_name, RNA_info, ATAC_info, ATAC_KS)
double_up_info <- double_up_info[order(double_up_info[,4]),]
saveRDS(double_up_info, file="double_down_2/double_up_info.rds")


################
# venn diagram #
################
library(VennDiagram)
tmp <- venn.diagram( list(RNA=unique(rna_up), ATAC=unique(peak_up)),
                     filename = NULL, height = 800, width = 800,  fill=c("yellow","green"), cex=0.5, cat.cex=0.5)
pdf("double_down_2/ovlp_venn.pdf")
grid.draw(tmp)
dev.off()


###############
# scatterplot #
###############
ovlp <- intersect(rownames(rna), names(median_ATAC_log2FC))
x <- rna[ovlp, "log2FC"]
y <- median_ATAC_log2FC[ovlp]
names(x) <- names(y) <- ovlp
pdf("double_down_2/scatter.pdf", 8, 8)
plot(x, y, pch=20, col="gray", xlim=c(-3,3), ylim=c(-2,2),
     xlab="RNA_log2FC (S/A)", ylab="ATAC_median_log2FC (S/A)")
cor(x, y)
points(x[double_up], y[double_up], col="red")
points(x[double_down], y[double_down], col="blue")
dev.off()

#######
# CDF #
#######
pdf("double_down_2/accessible_vs_expression.pdf")
plot(ecdf(rna$log2FC), pch=20, col="black", xlim=c(-2,2), main="expression change for differentially accessible genes", xlab="RNA log2FC (A/S)", ylab="CDF")
lines(ecdf(rna[peak_up,"log2FC"]), col="red", pch=20, cex=0.1)
lines(ecdf(rna[peak_down,"log2FC"]), col="blue", pch=20, cex=0.1)
abline(v=0, lty="dashed")
legend("topleft", legend = c("all genes","genes more accessible in adherent",
                             "genes less accessible in adherent"), fill=c("black","red","blue"))
p1 <- ks.test(rna[peak_up,"log2FC"], rna$log2FC, alternative="less")$p.value
p2 <- ks.test(rna[peak_down,"log2FC"], rna$log2FC, alternative="greater")$p.value
text(1, 0.1, paste0("accessibility+ vs. all p < ", signif(p1,3)))
text(1, 0.2, paste0("accessibility- vs. all p < ", signif(p2,3)))
dev.off()

#########
# genes #
#########
double_down_genes <- double_up
map <- select(org.Hs.eg.db, keys=double_down, columns="GENENAME", keytype="SYMBOL")
rownames(map) <- map[,1]
double_down_genes$full_name <- map[double_down_genes$gene,2]

double_down_atac <- do.call(rbind, lapply(1:length(double_down), function(i) {
  tmp_pks <- peaks[peaks$gene_name==double_down[i],]
  tmp_pks <- tmp_pks[which.min(tmp_pks$log2FC),]
}))
double_down_atac$KS_pvalue <- res[double_down_atac$gene_name,2]
gene_list <- data.frame(double_down_genes[,c("gene","full_name","log2FC","padj")], double_down_atac[,c("name","region","annot","distToNearest","KS_pvalue","log2FC","padj")])
write.xlsx(gene_list,"double_down_2/double_down_genes.xlsx", sheetName="all_genes")

###########
# heatmap #
###########
library(ComplexHeatmap)
library(circlize)
toplot <- counts_norm[double_up,c(4:6,1:3)]
labels <-c(rep("S",3),rep("A",3))
ha=HeatmapAnnotation(df = data.frame(type=as.character(labels)), 
                     col = list(type = c("S"="#cbd5e8","A" = "#fdcdac")))
pdf("double_down_2/Heatmap_double_down.pdf", 6, 8)
Heatmap(t(scale(t(log(toplot)))), show_row_names=T, top_annotation=ha, 
        show_row_dend=F, cluster_columns=F,  row_names_gp = gpar(fontsize = 10), 
        cluster_rows=T ,column_names_max_height=unit(6, "cm"))
dev.off()

###############################
# check for RUNX1 association #
###############################
# set all RUNX1 associated peak to log2FC=0
m <- readRDS("../motif_analysis/fimo/motif_hits.rds")
tmp <- atlas
tmp$RUNX1 <- 0
tmp[m[,"RUNX1"]==1]$RUNX1 <- 1
peaks <- data.frame(mcols(tmp))
peaks$name <- names(tmp)
peaks$gene <- peaks$gene_name
peaks$gene.dist <- atlas$distToNearest
peaks <- peaks[peaks$distToNearest < 50000, ]
genes <- read.xlsx2("double_down_2/double_down_genes.xlsx", sheetIndex=1, startRow=2)
peaks <- peaks[peaks$gene %in% as.character(genes$gene), ]
RUNX1_neg <- peaks[peaks$RUNX1==0, "log2FC"]
RUNX1_pos <- peaks[peaks$RUNX1==1, "log2FC"]

# expected # of RUNX1 peaks per gene
p <- sum(m[,"RUNX1"]==1) / nrow(m)
RUNX1_bound <- do.call(rbind, tapply(peaks$RUNX1, peaks$gene_name, function(x) {
  tmp <- c("pos"=sum(x==1), "neg"=sum(x==0))
  res <- binom.test(tmp, p=p, alternative="greater")$p.value
  return(c(tmp, binomial.p=res))
}))
RUNX1_bound <- cbind(RUNX1_bound, OR=RUNX1_bound[,1]/RUNX1_bound[,2]/p)
RUNX1_bound <- RUNX1_bound[order(RUNX1_bound[,4], decreasing=T),]
write.xlsx(RUNX1_bound, "double_down_2/RUNX1_bound.xlsx")


pdf("double_down_2/RUNX1_with_candidates.pdf")
plot(ecdf(atlas$log2FC), pch=20, col="gray", xlim=c(-2,2), main="RUNX1+ vs. RUNX1- peaks", xlab="accessibility log2FC (S/A)", ylab="CDF")
lines(ecdf(RUNX1_pos), col="blue", pch=20)
lines(ecdf(RUNX1_neg), col="cyan", pch=20)
abline(v=0, lty="dashed")
legend("topleft", legend = c("all pks","RUNX1- pks in candidates","RUNX1+ pks in candidates"), fill=c("gray","cyan","blue"))
p1 <- ks.test(RUNX1_neg, atlas$log2FC, alternative="greater")$p.value
p2 <- ks.test(RUNX1_pos, RUNX1_neg, alternative="greater")$p.value
text(1, 0.1, paste0("RUNX1- vs. all p < ", signif(p1,3)))
text(1, 0.2, paste0("RUNX1+ vs. RUNX1- p < ", signif(p2,3)))
dev.off()


###################################
# motif enrichment in double down #
###################################
m <- readRDS("../motif_analysis/fimo/motif_hits.rds")

# how does promoter peaks change?
# how does non-promoter peaks change?
pdf("double_down_2/promoter_vs_other_pks.pdf")
plot(ecdf(peaks[peaks$annot=="promoter","log2FC"]), col="green", pch=20, xlim=c(-2,2),
     xlab="accessibility log2FC (S/A)", main="accessibility of promoter vs. other pks in candidates")
lines(ecdf(peaks[peaks$annot!="promoter","log2FC"]), col="blue", pch=20)
abline(v=0, lty="dashed")
legend("topleft", legend = c("promoter","other"), fill=c("green","blue"))
dev.off()
# promoter peaks change less than enhancer peaks

# 1) motif enrichment in all peaks associated with double down genes

# KS test for each TF
tmp <- gsub("peak","Unit",peaks$name)
m_tmp <- m[tmp,]
res <- matrix(0, ncol(m), 3, dimnames=list(colnames(m), c("hits","stat","pvalue")))
for (i in 1:ncol(m)) {
  x = peaks[m_tmp[,i]==1,]$log2FC
  y = peaks[m_tmp[,i]==0,]$log2FC
  res[i,1] <- length(x)
  res[i,2] <- ks.test(x, y, alternative="greater")$statistic
  res[i,3] <- ks.test(x, y, alternative="greater")$p.value
  print(i)
}
res <- res[order(res[,3]), ]
write.xlsx2(res, file="double_down_2/motif_enrich_in_double_down.xlsx", sheetName="KS_test")

# odds ratio for each TF (hit/nonhit ratio within candidates, and hit/nohit ratio outside of candidates)
tmp <- as.numeric(gsub("peak","",peaks$name))
res <- matrix(0, ncol(m),5,dimnames=list(colnames(m), c("double_down_1","double_down_0","other_1","other_0","Odds_Ratio")))
for (i in 1:ncol(m)){
  a <- m[tmp,i]
  b <- m[-tmp, i]
  or <- (sum(a == 1) / sum(a == 0))/(sum(b == 1) / sum(b == 0))
  res[i,1] <- sum(a==1)
  res[i,2] <- sum(a==0)
  res[i,3] <- sum(b==1)
  res[i,4] <- sum(b==0)
  res[i,5] <- or
  print(i)
}
res <- res[order(res[,5], decreasing=T),]
write.xlsx2(res, file="double_down_2/motif_enrich_in_double_down.xlsx", sheetName="odds_ratio", append=T)

x <- res[,2]/(res[,2]+res[,3])
y <- res[,6]
pdf("double_down_2/top_enriched_TFs_in_42_genes.pdf")
plot(x,y, pch=20, xlab="% hits", ylab="odds ratio", main="top enriched TFs in 42 genes")
text(x[1:10], y[1:10], res[1:10,1], pos=1)
dev.off()

###########################
# RUNX1 targets vs. other #
###########################
runx1_bound <- unique(atlas[m[,"RUNX1"]==1]$gene_name)
runx1_log2FC <- rna[runx1_bound,"log2FC"]
other_log2FC <- rna[!rownames(rna)%in%runx1_bound,"log2FC"]

pdf("double_down_2/RUNX1_target_expression_change.pdf")
plot(ecdf(runx1_log2FC), pch=20, col="blue", xlim=c(-1,1), main="RUNX1+ genes sv. RUNX1- genes", xlab="expression log2FC (S/A)", ylab="CDF")
lines(ecdf(other_log2FC), pch=20, col="gray")
pval <- ks.test(runx1_log2FC, other_log2FC, alternative = "greater")$p.value
text(0.5, 0.2, paste0("KS test p < ", signif(pval, 3)))
dev.off()


######################
# 42 genes as a plot #
######################
# this doesn't work
a <- read.xlsx2("double_down_2/RUNX1_bound_rev.xlsx", sheetIndex=1)
a[,2] <- as.numeric(as.character(a[,2]))
a[,3] <- as.numeric(as.character(a[,3]))
x <- a[,2]
y <- as.numeric(as.character(a[,5]))
pdf("double_down_2/42_candidates_by_RUNX1.pdf")
plot(x,y,pch=20,xlab="RUNX1 binding sites",ylab="odds ratio")
text(x,y,as.character(a[,1]),pos=2)
dev.off()

