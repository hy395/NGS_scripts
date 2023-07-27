## trim and align SE reads with bowtie2
## this is for RNA-seq data, so I set "remove duplicate" to "no"
## adapter sequence for trimming should be provided in pwd with the name "adapter.fa"
## this is updated on 11.20.2017

INDIR=$1
OUTDIR=$2


mkdir ${OUTDIR}
mkdir ${OUTDIR}/trim_FASTQC
mkdir star_alignment

for ENTRY in ${INDIR}/*.gz
do
 NAME=`basename ${ENTRY} .fastq.gz`
 fastqc -q -t 8 -o ${OUTDIR}/trim_FASTQC ${OUTDIR}/${NAME}_trim.fastq.gz
 /home/hy395/programs/myscripts/fastq2bam_SE_star.sh ${OUTDIR}/${NAME}_trim.fastq.gz /data/leslie/genome_index/hg19.star/ star_alignment 8
 # for spikein
 #/ifs/e63data/leslielab/hy395/executable/myscripts/fastq2bam_SE_star.sh ${OUTDIR}/${NAME}_trim.fastq.gz /ifs/e63data/leslielab/hy395/hg19_spikein.star/genome star_alignment 10
 #rm ${OUTDIR}/${NAME}_trim.fastq.gz
 echo "$NAME"
done
