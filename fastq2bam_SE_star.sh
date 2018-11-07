#! /usr/bin/env bash
## argments:
## $1: fastq file 
## $2: sample name
## $3: output dir
## $4: adapter
## $5: star index
## $6: dir for star output stats
## $7: # threads

## this is aligning SE reads with STAR (RNA-seq only). resulting bam files are filtered, but NOT duplicate-removed. updated on 4.18.2017
## outputs are in the same folder as inputs.
## filtering: -F 0x204 removes 1) reads that don't map. 2) QC failed reads.

FILE=$1
SAMPLE=$2
OUTDIR=$3
ADAPTER=$4
GENOME_DIR=$5
STARDIR=$6
CORES=$7

# trim
java -jar /home/hy395/programs/source/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads ${CORES} -phred33 ${FILE} ${OUTDIR}/${SAMPLE}_trim.fastq.gz ILLUMINACLIP:${ADAPTER}:2:30:7 SLIDINGWINDOW:3:15 MINLEN:36

# alignment
STAR --genomeLoad NoSharedMemory --genomeDir ${GENOME_DIR} \
--readFilesIn ${OUTDIR}/${SAMPLE}_trim.fastq.gz \
--readFilesCommand zcat \
--runThreadN ${CORES} \
--limitBAMsortRAM 20000000000 \
--outFileNamePrefix ${STARDIR}/${SAMPLE} \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNmax 5 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outStd BAM_SortedByCoordinate \
--outSAMstrandField intronMotif > ${OUTDIR}/${SAMPLE}.bam

# filtering
BAM=${OUTDIR}/${SAMPLE}.bam
samtools view -h ${BAM}|samtools view -bh -F 0x204 -q 10 -S - > ${BAM/.bam/.filter.bam}
samtools index ${BAM/.bam/.filter.bam}
samtools flagstat ${BAM/.bam/.filter.bam} > ${BAM/.bam/.filter.flagstat}

rm ${BAM}
rm ${OUTDIR}/${SAMPLE}_trim.fastq.gz
