#############
# download MDS microarray
#############
library(GEOquery)
library(data.table)
library(preprocessCore)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(hgu133a.db)

##############################
# get microarray data matrix #
##############################

# use getGEO
# important to do quantile normalization
GSE <- getGEO("GSE19429",GSEMatrix=TRUE)
exp <-GSE[[1]]
exp.m <- exprs(exp)
temp <- normalize.quantiles(exp.m)
dimnames(temp) <- dimnames(exp.m)
exp.m <- temp
probe2symbol <- select(hgu133plus2.db, keys= rownames(exp.m), column = c("SYMBOL"), keytype = "PROBEID")
probe2symbol <- probe2symbol[!is.na(probe2symbol[,2]),] # remove probes with no symbol match
probe2symbol <- probe2symbol[!duplicated(probe2symbol[,1]),]
rownames(probe2symbol) <- probe2symbol[,1]
exp.m <- exp.m[rownames(exp.m)%in%probe2symbol[,1],] # filter expression matrix
# create a data table and summary to gene level signal (median)
exp.dt <- data.table(exp.m, symbol=probe2symbol[rownames(exp.m),2])
exp.dt <- exp.dt[,lapply(.SD,median),by=symbol]
names <- exp.dt$symbol
exp.dt[,symbol:=NULL]
exp.m <- data.frame(exp.dt)
rownames(exp.m) <- names
# phenotype
pheno <- pData(phenoData(exp))
saveRDS(list(x=exp.m,y=pheno),file="data/mds_microarray.rds")

# get VLK microarray
GSE <- getGEO("GSE1159",GSEMatrix=TRUE,getGPL=F)
exp <-GSE[[1]]
exp.m <- exprs(exp)
temp <- normalize.quantiles(exp.m)
dimnames(temp) <- dimnames(exp.m)
exp.m <- temp
probe2symbol <- select(hgu133a.db, keys= rownames(exp.m), column = c("SYMBOL"), keytype = "PROBEID")
probe2symbol <- probe2symbol[!is.na(probe2symbol[,2]),] # remove probes with no symbol match
probe2symbol <- probe2symbol[!duplicated(probe2symbol[,1]),]
rownames(probe2symbol) <- probe2symbol[,1]
exp.m <- exp.m[rownames(exp.m)%in%probe2symbol[,1],] # filter expression matrix
# create a data table and summary to gene level signal (median)
exp.dt <- data.table(exp.m, symbol=probe2symbol[rownames(exp.m),2])
exp.dt <- exp.dt[,lapply(.SD,median),by=symbol]
names <- exp.dt$symbol
exp.dt[,symbol:=NULL]
exp.m <- data.frame(exp.dt)
rownames(exp.m) <- names
# phenotype
pheno <- pData(phenoData(exp))
pheno <- pheno[,1:2]
pheno[,1] <- gsub("#","",as.character(pheno[,1]))
colnames(exp.m) <- pheno[colnames(exp.m),1]
saveRDS(exp.m,file="data/vlk.rds")

###############
# download TCGA AML data
###############
library(org.Hs.eg.db)
library(data.table)

files <- list.files(pattern="*.counts",path="data/tcga_aml",recursive=T, full.names=T)
counts.l <- lapply(files, function(x) {
  temp <- read.table(x, colClasses=c("character","numeric"))
  return(temp[1:60483,]) # last few lines are outputs from ht-seq
})
names <- counts.l[[1]][,1]
counts.l <- lapply(counts.l, function(x) x[,2])
counts <- do.call(cbind, counts.l)
rownames(counts) <- names
tomap <- gsub("\\..*","",rownames(counts))
ensembl2symbol <- select(org.Hs.eg.db, keys= tomap, column = c("SYMBOL"), keytype = "ENSEMBL")
ensembl2symbol <- ensembl2symbol[!is.na(ensembl2symbol[,2]),] # remove probes with no symbol match
ensembl2symbol <- ensembl2symbol[!duplicated(ensembl2symbol[,1]),]
rownames(ensembl2symbol) <- ensembl2symbol[,1]
counts <- counts[tomap%in%ensembl2symbol[,1],] # filter count matrix
rownames(counts)  <- ensembl2symbol[gsub("\\..*","",rownames(counts)),2]
colnames(counts) <- paste0("tcga_aml_",1:ncol(counts))
saveRDS(counts, file="data/tcga_aml_rnaseq.rds")


####################
# download TCGA AML data  from firehose
####################
fh_aml <- read.table("data/LAML.rnaseq__illuminaga_rnaseq__bcgsc_ca__Level_3__gene_expression__data.data.txt",sep="\t", colClasses="character")
fh_aml_count <- fh_aml[,c(1,grep("raw",fh_aml[2,]))]
colnames(fh_aml_count) <- fh_aml_count[1,]
fh_aml_count <- fh_aml_count[-c(1:2,grep("\\?",fh_aml_count[,1])),]
fh_aml_count[,1] <- gsub("\\|.*","",fh_aml_count[,1])
genes <- fh_aml_count[,1]
fh_aml_count <- apply(fh_aml_count[,-1], 2, function(x) {
  return(as.numeric(x))
  })
count.dt <- data.table(genes=genes,fh_aml_count)
count.dt <- count.dt[,lapply(.SD,sum),by=genes]
names <- count.dt$genes
count.dt[,genes:=NULL]
fh_aml.m <- data.frame(count.dt)
rownames(fh_aml.m) <- names
saveRDS(fh_aml.m, file="data/firehose_aml.rds")



#####################################
# count read for normal bone marrow #
#####################################
library(rtracklayer)
library(DESeq2)
library(GenomicAlignments)

#hg19
annotations <- import("../genome/ENSEMBL.homo_sapiens.release-75/Homo_sapiens.GRCh37.75.gtf", format="gtf")
saveRDS(annotations,file="../genomes/Homo_sapiens.GRCh37.75.rds")

chrs<-c(1:22,"X","Y")
#annotation GRList
annotations<-readRDS("../genomes/Homo_sapiens.GRCh37.75.rds")
annotations<-keepSeqlevels(annotations,chrs)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
genes <- unique(annotations$gene_name)
exons <- annotations[annotations$type == 'exon']
exons.genes <- disjoin(split(exons, factor(exons$gene_name,levels=genes)))

#rRNA GRList
tRNA_rRNA <- import("../genomes/hg19_trna_rrna.gtf", format="gtf")
tRNA_rRNA<-keepSeqlevels(tRNA_rRNA,paste0("chr",chrs))
rrnas <-unique(tRNA_rRNA$gene_id)
tRNA_rRNA <- disjoin(split(tRNA_rRNA, factor(tRNA_rRNA$gene_id,levels=rrnas)))
tRNA_rRNA <- endoapply(tRNA_rRNA, function (x) {mcols(x) <- NULL; x})

#removing overlapping regions using find overlap
overlap<-findOverlaps(exons.genes,tRNA_rRNA)
rrnas<-unique(overlap@queryHits)
exons.genes<-exons.genes[-rrnas]
bams <- dir("data/external/nbm/", pattern = "*.bam$", recursive = TRUE,
            include.dirs = TRUE, full.name = TRUE)
counts <- summarizeOverlaps(exons.genes, bams, mode = "IntersectionNotEmpty",
                            ignore.strand = TRUE, singleEnd=TRUE)
saveRDS(counts,file="data/counts_nbm.rds")


###################
# IPSC fibroblast #
###################
ipsc_fb <- lapply(list.files("iPSC_arm_fibroblast/","*.tsv",full.names=T), function(x) {
  temp <- fread(x)
  cnts <- temp$posterior_mean_count
  names(cnts) <- temp$gene_id
  return(cnts)  
})
ipsc_fb <- Reduce(cbind, ipsc_fb)
rownames(ipsc_fb)
symbols<-select(org.Hs.eg.db, keys=gsub("\\..*","",rownames(ipsc_fb)),columns= c("SYMBOL","ENSEMBL"), keytype="ENSEMBL")
symbols <- symbols[!duplicated(symbols[,1]),]
symbols <- symbols[!is.na(symbols[,2]),]
rownames(symbols) <- symbols[,1]
rownames(ipsc_fb) <- gsub("\\..*","",rownames(ipsc_fb))
ipsc_fb <- ipsc_fb[rownames(symbols), ]
rownames(ipsc_fb) <- symbols[rownames(ipsc_fb),2]
colnames(ipsc_fb) <- c("ipscfb_1","ipscfb_2")





