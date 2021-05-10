#! /bin/bash
## bowtie2 alignment for PE reads
## duplicate reads are removed by picard
## adapter file has to be supplied

## softwares:
# cutadapt
# bowtie2
# samtools
# picard
# use 8G of memory

## 1 is Bowtie2 index file
## 2 is R1.fastq.gz, before trimming
## 3 is R2.fastq.gz, before trimming
## 4 is the adapter file
## 5 sample name, can include output directory
## 6 is number of cores

INDEX=$1
R1=$2
R2=$3
ADAPTER=$4
SAMPLE=$5
CORES=$6

# trim: need to take adapter. use cutadapt instead of trimmomatic. trimmomatic throw away entire read when detected whole adapter. (reported in v0.36)
cutadapt -a file:${ADAPTER} -A file:${ADAPTER} -o ${SAMPLE}.R1.fq -p ${SAMPLE}.R2.fq --pair-filter=any --minimum-length=30 ${R1} ${R2}
#java -jar /home/yuanh/programs/source/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog ${SAMPLE}.trimlog ${R2} ${R2} ${SAMPLE}_R1_paired.fq.gz ${SAMPLE}_R1_unpaired.fq.gz ${SAMPLE}_R2_paired.fq.gz ${SAMPLE}_R2_unpaired.fq.gz ILLUMINACLIP:${ADAPTER}:2:30:10:2:keepBothReads MINLEN:36
#zcat ${SAMPLE}_R1_paired.fq.gz > ${SAMPLE}.R1.fq & zcat ${SAMPLE}_R2_paired.fq.gz > ${SAMPLE}.R2.fq

# alignment: dont align discordant reads, dont align single reads. For ATAC-seq, need to change the default -X(--maxins) option to 1000. This is the maximum frag length. The default is 500.
# filtering: -f 0x02 only keep reads mapped in proper pair. -q 20 only keep fragments where the total q of two mates are > 20.
bowtie2 -x ${INDEX} -1 ${SAMPLE}.R1.fq -2 ${SAMPLE}.R2.fq --no-mixed --no-discordant -p ${CORES} -X 1000 | samtools view -S -bf 0x02 -q 20 -o ${SAMPLE}.nonSorted.bam -
rm ${SAMPLE}.R1.fq ${SAMPLE}.R2.fq

# sort the bam file
samtools sort -m 4G -o ${SAMPLE}.filter.bam ${SAMPLE}.nonSorted.bam

# remove duplicates
java -Xmx8g -jar /home/yuanh/programs/source/picard-2.20.3/picard.jar MarkDuplicates REMOVE_DUPLICATES=true AS=true I=${SAMPLE}.filter.bam O=${SAMPLE}.filter.rmdup.bam M=${SAMPLE}_rmdup_metric.txt

# index and clean up
samtools index ${SAMPLE}.filter.rmdup.bam
samtools flagstat ${SAMPLE}.filter.rmdup.bam > ${SAMPLE}.flagstat
rm ${SAMPLE}.nonSorted.bam ${SAMPLE}.filter.bam
