library(SRAdb)
sqlfile <- '~/programs/genomes/SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
conversion <- sraConvert( 'SRP049300', sra_con = sra_con )
dbListFields(sra_con,"experiment") # "experiment" has description for each EXP
rs <- dbGetQuery(sra_con, paste( "select experiment_ID, title, experiment_alias, experiment_accession from experiment where",
                                 "submission_accession like 'SRA196613'",sep=" "))
rownames(rs) <- rs$experiment_accession
conversion$description <- rs[conversion$experiment,]$title

# decide which to download
conversion <- conversion[order(conversion$run),]
todownload <- conversion

# annotations of downloaded files
dir.create("raw")
write.csv(todownload, file="raw/files.txt")

# download through curl
todownload <- listSRAfile( todownload$run, sra_con, fileType = 'sra' )
write.table(todownload$ftp, file="urls.txt", col.names = F, row.names = F, quote = F)
setwd("raw")
system("xargs -n 1 curl -O -L <../urls.txt")

# convert them to fastq
system("for i in *.sra; do fastq-dump --gzip ${i}; done")
system("rm *.sra")

