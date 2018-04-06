## trim and align SE reads with bowtie2
## adapter sequence for trimming should be provided in pwd with the name "adapter.fa"
## this is updated on 9.20.2017
## adapted to lilac


INDIR=$1
OUTDIR=$2
g_index=$3

mkdir ${OUTDIR}
mkdir ${OUTDIR}/trim_FASTQC

for ENTRY in ${INDIR}/*.gz
do
 NAME=`basename ${ENTRY} .fastq.gz`
 java -jar /home/hy395/programs/source/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 ${ENTRY} ${OUTDIR}/${NAME}_trim.fastq.gz ILLUMINACLIP:adapter.fa:2:30:7 SLIDINGWINDOW:3:15 MINLEN:36
 fastqc -q -t 8 -o ${OUTDIR}/trim_FASTQC ${OUTDIR}/${NAME}_trim.fastq.gz
 /home/hy395/programs/myscripts/fastq2bam_SE_bowtie.sh ${OUTDIR}/${NAME}_trim.fastq.gz ${g_index} 8
 rm ${OUTDIR}/${NAME}_trim.fastq.gz
 echo "$NAME"
done
