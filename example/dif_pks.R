library(rtracklayer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(DESeq2)
source("../../BRD/functions.R")
dir.create("dif_pks")

# reading data
load("../safe_harbor/safe_harbor.rdt")
atlas <- readRDS("rdata/atac_atlas.rds")
names(atlas) <- atlas$name
atlas <- add_annot(atlas)

# read counts on atlas
inputPath <-list.files("rdata", pattern="rmdup", full.names=T)
bams <- lapply(inputPath, readRDS)
names(bams) <- gsub("rdata/","",gsub("_IGO_.*","",inputPath))
atlas <- countReads(atlas, bams)

# DESeq
cnts <- data.frame(mcols(atlas)[,6:12])
celltype <- gsub(".*_","",colnames(cnts))
sample_id <- gsub("_.*","",colnames(cnts))
colnames(cnts) <- paste0(celltype, "_", sample_id)
dds <- DESeqDataSetFromMatrix(cnts, colData=data.frame(celltype=celltype,sample_id=sample_id), design=~celltype)
dds <- estimateSizeFactors(dds)
model <- DESeq(dds, betaPrior = F)
res <- results(model, alpha = 0.05, contrast=c("celltype","CD4","CD8"), independentFiltering = F)
atlas$log2FC <- res$log2FoldChange
atlas$padj <- res$padj

# differential pks
dif_pks  <- atlas[atlas$padj < 0.05]
dif_pks <- dif_pks[order(dif_pks$log2FC)]
dif_pks$region <- paste0(seqnames(dif_pks),":",start(dif_pks),"-",end(dif_pks))
igv_window <- resize(dif_pks, 2000, fix="center")
dif_pks$igv_window <- paste0(seqnames(igv_window),":",start(igv_window),"-",end(igv_window))
write.xlsx(data.frame(mcols(dif_pks)), file="dif_pks/dif_pks.xlsx", sheetName="sheet1")

