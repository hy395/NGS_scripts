# DESeq2 functions
# 5.24.2017 remove the RPKM function, use DESeq2::fpkm instead
library(org.Hs.eg.db)
library(GOstats)

#############
# symbol2id #
#############
symb2id <- function(temp) {
  if (length(temp) == 0) {
    return(NULL)
  } else {
    a <- select(org.Hs.eg.db, keys=temp,columns= c("ENTREZID"),keytype="SYMBOL")$ENTREZID
    return(unique(a[!is.na(a)]))
  }
}

########################
# filter based on fpkm #
########################
fpkmFilterRange<-function(dds,minSample,fpkmRange){
  #rpkmRange: consider minRPKM chosen from this range.
  fpkm<-assays(dds)$fpkm
  genes.left<-vector()
  for (i in fpkmRange){
    rows<-rowSums(fpkm>=i)>=minSample
    genes.left<-c(genes.left,sum(rows))
  }
  return( cbind(fpkmRange,genes.left) )
}
fpkmFilter<-function(dds,minSample,minFPKM){
  #minRPKM and minSample: consider genes that have RPKM>=minRPKM in at least minSample samples.
  fpkm<-assays(dds)$fpkm
  rows<-rowSums(fpkm>=minFPKM)>=minSample
  cat(sprintf("total genes: %i\nremaining genes: %i\n",length(rows),sum(rows)))
  return(dds[rows,])
}

####################################
# filter based on normalized count #
####################################
countFilterRange<-function(DESeqDataSet,minSample,countRange){
  counts<-counts(DESeqDataSet,normalized=T)
  genes.left<-vector()
  for (i in countRange){
    rows<-rowSums(counts>=i)>=minSample  
    genes.left<-c(genes.left,sum(rows))
  }
  return( cbind(countRange,genes.left) )
}
countFilter<-function(DESeqDataSet,minSample,minCount){
  counts<-counts(DESeqDataSet,normalized=T)
  rows<-rowSums(counts>=minCount)>=minSample
  cat(sprintf("total genes: %i\nremaining genes: %i\n",length(rows),sum(rows)))
  return (DESeqDataSet[rows,] )
}

################
# volcano plot #
################
volcano<-function(res.df,log2FoldThres,padjThres,sig.genes=NULL,main){
  #sig.genes is the genes that I want to label, in ensembl_gene_id
  x<-res.df$log2FoldChange
  y<-res.df$padj
  plot(x,-log10(y),main=main,xlab="",ylab="-log10 Padj",pch=20)
  downRows<-x<(-log2FoldThres)&y<padjThres;downRows[is.na(downRows)]<-FALSE
  upRows<-x>log2FoldThres&y<padjThres;upRows[is.na(upRows)]<-FALSE
  mtext(paste("total genes: ",length(downRows),"; up-regulated genes: ",sum(upRows,na.rm=T),"; down-regulated genes: ",sum(downRows,na.rm=T),sep=""), side=1,line=2)
  cat(sprintf("total genes: %i\nup-regulated genes: %i\ndown-regulated genes: %i\n",length(downRows),sum(upRows,na.rm=T),sum(downRows,na.rm=T)))
  points(x[downRows],-log10(y[downRows]),col="blue",pch=20)
  points(x[upRows],-log10(y[upRows]),col="red",pch=20)
  #if fed sig.genes list
  if (!is.null(sig.genes)){
    #plot signature genes
    sigRows<-rownames(res.df)%in%sig.genes&(x<(-log2FoldThres)|x>log2FoldThres)&y<padjThres
    if (sum(sigRows,na.rm=T)>0){
      text(x[sigRows],-log10(y[sigRows]),res.df[sigRows,7],pos=3,cex=1.2)
    }
    downSigRows<-rownames(res.df)%in%sig.genes&x<(-log2FoldThres)&y<padjThres;downSigRows[is.na(downSigRows)]<-FALSE
    upSigRows<-rownames(res.df)%in%sig.genes&x>log2FoldThres&y<padjThres;upSigRows[is.na(upSigRows)]<-FALSE
    mtext(paste("total signature genes: ",sum(rownames(res.df)%in%sig.genes),"; up-regulated genes: ",sum(upSigRows,na.rm=T),"; down-regulated genes: ",sum(downSigRows,na.rm=T),sep=""), side=1,line=3)
    cat(sprintf("total signature genes: %i\nup-regulated signature genes: %i\ndown-regulated signature genes: %i\n",sum(rownames(res.df)%in%sig.genes),sum(upSigRows,na.rm=T),sum(downSigRows,na.rm=T)))    
    down.sig.genes<-res.df[downSigRows,]
    up.sig.genes<-res.df[upSigRows,]
  }
  #output
  down.genes<-res.df[downRows,]
  up.genes<-res.df[upRows,]
  if (!is.null(sig.genes)){
    output<-list(down.genes=down.genes,up.genes=up.genes,down.sig.genes=down.sig.genes,up.sig.genes=up.sig.genes)    
  } else {
    output<-list(down.genes=down.genes,up.genes=up.genes)
  }
  return(output)
}

###############
# GO analysis #
###############
runGOstats <- function(genes, genes.all, conditional = TRUE, p.value =
                         0.05, ontology = c('BP', 'MF', 'CC'),
                       annotation = 'org.Hs.eg.db',
                       testDirection = 'over', KEGG=FALSE) {
  # Run a single GO enrichment test
  if (length(genes) == 0) {
    return (NULL)
  }
  if (all(is.na(select(org.Hs.eg.db, keys = genes, columns = "ONTOLOGY", keytype = "ENTREZID")$ONTOLOGY))) {
    return(NULL)
  }
  ontology <- match.arg(ontology)
  genes <- as.character(genes)
  genes.all <- as.character(genes.all)
  if (KEGG==TRUE){
    params <- new("KEGGHyperGParams", geneIds = genes,
                  universeGeneIds = genes.all, annotation = annotation,
                  pvalueCutoff = p.value ,testDirection = testDirection)
    go.res <- hyperGTest(params)
    go.res
  } else {
    params <- new("GOHyperGParams", geneIds = genes,
                  universeGeneIds = genes.all, annotation = annotation,
                  ontology = ontology, pvalueCutoff = p.value,
                  conditional = conditional, testDirection =
                    testDirection)
    go.res <- hyperGTest(params)
    go.res
  }
}

#create my summary of GO result
mySummary<-function(GOresult,minCount){
  if (is.null(GOresult)) {
    return(NA)
  }
  ids.symbols<-select(org.Hs.eg.db, keys=GOresult@geneIds,columns= c("SYMBOL"),keytype="ENTREZID")
  GOresult.summ<-summary(GOresult)
  gos<-GOresult.summ[,1]
  a<-select(org.Hs.eg.db, keys=gos,columns= c("SYMBOL"),keytype="GOALL")
  a.list<-split(a,a$GOALL)
  genes<-vector()
  for (i in 1:length(a.list)){
    overlap<-intersect(unique(a.list[[i]]$SYMBOL),ids.symbols[,2])
    overlap.collapse<-paste(overlap,collapse=",")
    genes[i]<-overlap.collapse
  }
  genes.collapse<-cbind(names(a.list),genes)
  rownames(GOresult.summ)<-GOresult.summ$GOBPID
  GOresult.summ[genes.collapse[,1],"genes"]<-genes.collapse[,2]
  GOresult.summ[GOresult.summ$Count>=minCount,]
}

#MAplot, volcano plot, list of peaks go analysis
deseq_analysis<-function(tf,subdds,factors,up,down,log2fThres=2,padjThres=1,GO=T){
  colData(subdds)$condition<-factor(factors)
  subdds <- DESeq(subdds, betaPrior=FALSE)
  result<-results(subdds,independentFiltering=FALSE,contrast=c("condition",up,down))
  #MAplot
  comparison<-paste0(down,"vs",up)
  png(paste0("figures/MAplot_",comparison,".png"));plotMA(result,alpha=0.05,main=comparison);dev.off()
  #Volcano plot
  png(paste0("figures/Volcano_",comparison,".png"))
  vol<-volcano(data.frame(result),log2FoldThres=log2fThres,padjThres=padjThres,main=comparison)
  dev.off()
  #get gene symbols of differentially expressed ones
  all_symbols <- combined.peaks[rownames(result)]$symbol
  up_symbols <- combined.peaks[rownames(vol[["up.genes"]])]$symbol
  down_symbols<-combined.peaks[rownames(vol[["down.genes"]])]$symbol
  #GO analysis
  #change gene symbols to entrez ids
  if(GO){
    all.id<-select(org.Hs.eg.db, keys=all_symbols,columns="ENTREZID",keytype="SYMBOL");all.id<-as.numeric(all.id[,2])
    up.table <- "empty"
    if(length(up_symbols)!=0){
      up.id<-select(org.Hs.eg.db, keys=up_symbols,columns="ENTREZID",keytype="SYMBOL");up.id<-as.numeric(up.id[,2])      
      if(length(up.id)>5){
        up.result<-runGOstats(up.id,all.id)
        up.table<-mySummary(up.result,2)
      }      
    }
    down.table <- "empty"
    if(length(down_symbols)!=0){
      down.id<-select(org.Hs.eg.db, keys=down_symbols,columns="ENTREZID",keytype="SYMBOL");down.id<-as.numeric(down.id[,2])            
      if(length(down.id)>5){
        down.result<-runGOstats(down.id,all.id)
        down.table<-mySummary(down.result,2)    
      }
    }
    filename<-paste0(tf,"_GO.xlsx")
    write.xlsx(data.frame(up.table),file=filename,sheetName=paste0(comparison,"_up"),append=file.exists(filename))
    write.xlsx(data.frame(down.table),file=filename,sheetName=paste0(comparison,"_down"),append=file.exists(filename))    
  }
  return(list(result=result,vol=vol,down_symbols=down_symbols,up_symbols=up_symbols))
}

# BROAD GSEA
# function for running GSEA
genGseaFiles <- function(dds_tmp, labels, outdir){
  # labels have to be factors
  # dir need to have / at the end of it.
  dir.create(outdir)
  dds_tmp <- dds_tmp[,order(labels)]
  labels <- sort(labels)
  line1 <- paste(c(length(labels),length(levels(labels)),1),collapse=" ")
  line2 <- paste(c("#",levels(labels)), collapse=" ")
  line3 <- paste(labels,collapse=" ")
  file.create(paste0(outdir,"labels.cls"), recursive=T)
  fileConn <- file(paste0(outdir,"labels.cls"),"w")
  writeLines(line1,fileConn)
  writeLines(line2,fileConn)
  writeLines(line3,fileConn)
  close(fileConn)
  # produce matrices for GSEA
  cnt_norm <- counts(dds_tmp,normalized=T)
  exprs <- data.frame(NAME=rownames(cnt_norm), description="na",cnt_norm)
  write.table(exprs, file=paste0(outdir,"expr.txt"), sep="\t", col.names=T, row.names=F, quote=F)
}

# 
runGSEA <- function(label_file, expr_file, gmt_file, compare, mem_size, permute_type, rpt_label, out_name){
  # label_file: .cls file.
  # expr_file: expression file.
  # gmt_file: gene set file.
  # compare: comparison. #A_versus_B
  # mem_size: memory size (in MB)
  # permute_type: phenotype or gene_set.
  #  Phenotype. This is the recommended method when there are at least seven (7) samples in each phenotype. # 25% FDR cutoff for significance.
  #  Gene_set. When you have fewer than seven (7) samples in any phenotype. # 5% FDR cutoff for signficance.
  # rpt_label: name of report.
  # out_name: name of output dir.
  cmd <- paste0("java -cp /home/hy395/programs/bin/gsea-3.0.jar -Xmx" , mem_size, "m xtools.gsea.Gsea -res ",expr_file,
                " -cls ", label_file, compare,
                " -gmx ",gmt_file,
                " -collapse false -mode Max_probe -norm None -nperm 1000 -permute ", permute_type,
                " -rnd_type no_balance -scoring_scheme weighted -rpt_label ", rpt_label,
                " -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ",out_name,
                " -gui false")
  system(cmd)  
}


###########
# Heatmap #
###########
countBin<-function(bam,peakBin, window_size){
  vectorized_m<-countOverlaps(peakBin,bam)
  output<-matrix(vectorized_m,ncol=window_size,byrow=T)
  return(output)
}
generateMatrix <- function(regions, beds, covs, window_size=1000) {
  peaks_centered_tiles_gr<- unlist(tile(regions,window_size))
  
  # rpm count for each bin
  bin_matrix <- mclapply(beds, function(x) countBin(x,peaks_centered_tiles_gr, window_size=window_size), mc.cores=40)
  for (i in 1:length(bin_matrix)){
    bin_matrix[[i]]<-bin_matrix[[i]]/covs[i]*1e6
  }
  # put a cap at the intensity
  #quantSig <- quantile(sapply(bin_matrix,function(x) apply(x,1,max)))
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

