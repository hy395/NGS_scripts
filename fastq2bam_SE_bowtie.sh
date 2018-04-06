#! /usr/bin/env bash
## Argments:
## $1: fastq file
## $2: bowtie2 index file with path
## $3: threads to be used in alignment

## this is aligning SE reads with bowtie2. resulting bam files are duplicate removed and filtered. updated on 4.18.2017.
## outputs are in the same folder as inputs.
## adopted from https://raw.githubusercontent.com/shenlab-sinai/chip-seq_preprocess/master/bin/fastq2bam_by_bowtie2.sh
## filtering: -F 0x204 removes 1) reads that don't map. 2) QC failed reads.

FILE=$1
BOWTIE_INDEX=$2
FILENAME=$(basename "$FILE")
FQDIR=$(dirname "$FILE")
EXT="${FILENAME##*.}"
FILENAME_BASE="${FILENAME%.*}"
SAM=${FQDIR}/${FILENAME_BASE}.sam
CORES=$3

case "$EXT" in
	fq | fastq | FQ | FASTQ ) bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} ${FILE} -S ${SAM}
	    ;;
	gz | GZ ) zcat ${FILE} | bowtie2 -p ${CORES} -x ${BOWTIE_INDEX} -S ${SAM} -
	    ;;
esac

samtools view -Sb ${SAM} > ${SAM/sam/nonSorted.bam}
samtools sort -m 5G -o ${SAM/sam/bam} ${SAM/sam/nonSorted.bam}
samtools view -h ${SAM/sam/bam}|grep -v "XS:"|samtools view -bh -F 0x204 -q 10 -S - > ${SAM/sam/filter.bam}
samtools rmdup -s ${SAM/sam/filter.bam} ${SAM/sam/filter.nodup.bam}
samtools index ${SAM/sam/filter.nodup.bam}
rm ${SAM} ${SAM/sam/nonSorted.bam} ${SAM/sam/bam} ${SAM/sam/filter.bam} 

