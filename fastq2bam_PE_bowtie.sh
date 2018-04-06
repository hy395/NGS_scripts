#! /bin/bash
## bowtie2 alignment for PE reads
## duplicate reads are removed by picard
## adapter file has to be supplied

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

# trim: need to take adapter.
java -jar /home/hy395/programs/source/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${R1} ${R2} ${SAMPLE}_R1_paired.fq.gz ${SAMPLE}_R1_unpaired.fq.gz ${SAMPLE}_R2_paired.fq.gz ${SAMPLE}_R2_unpaired.fq.gz ILLUMINACLIP:${ADAPTER}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# alignment: dont align discordant reads, dont align single reads. For ATAC-seq, need to change the default -X(--maxins) option to 1000. This is the maximum frag length. The default is 500.
# filtering: -f 0x02 only keep reads mapped in proper pair. -q 20 only keep fragments where the total q of two mates are > 20.
zcat ${SAMPLE}_R1_paired.fq.gz > ${SAMPLE}.R1.fq & zcat ${SAMPLE}_R2_paired.fq.gz > ${SAMPLE}.R2.fq
bowtie2 -x ${INDEX} -1 ${SAMPLE}.R1.fq -2 ${SAMPLE}.R2.fq --no-mixed --no-discordant -p ${CORES} -X 1000 | samtools view -S -bf 0x02 -q 20 -o ${SAMPLE}.nonSorted.bam -
rm ${SAMPLE}.R1.fq ${SAMPLE}.R2.fq
rm ${SAMPLE}_R1_paired.fq.gz ${SAMPLE}_R1_unpaired.fq.gz ${SAMPLE}_R2_paired.fq.gz ${SAMPLE}_R2_unpaired.fq.gz

# sort the bam file
samtools sort -m 5G -o ${SAMPLE}.filter.bam ${SAMPLE}.nonSorted.bam

# remove duplicates
java -jar /home/hy395/programs/source/picard-tools-1.141/picard.jar MarkDuplicates REMOVE_DUPLICATES=true AS=true I=${SAMPLE}.filter.bam O=${SAMPLE}.filter.rmdup.bam M=${SAMPLE}_rmdup_metric.txt

# index and clean up
samtools index ${SAMPLE}.filter.rmdup.bam
samtools flagstat ${SAMPLE}.filter.rmdup.bam > ${SAMPLE}.flagstat
rm ${SAMPLE}.nonSorted.bam ${SAMPLE}.filter.bam

