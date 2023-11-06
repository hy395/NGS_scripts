#!/bin/bash

## this is aligning PE reads with STAR (RNA-seq only). 
# resulting bam files are filtered, but NOT duplicate-removed. updated on 7.26.23

#############
# arguments #
#############

# param defaults
ADAPTER='/home/yuanh/programs/source/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa' # default for Illumina RNA-seq
GENOME_DIR='/scratch4/yuanh/star_index/star41'
STAR='/home/yuanh/programs/source/STAR-2.7.3a/bin/Linux_x86_64_static/STAR'
OUTDIR='bams'
STAR_out='STAR_out'
CORES=16
ramGB=20

function show_usage() {
    echo "Usage: ./fastq2bam_PE_star.sh --r1 <_R1.fq.gz> --r2 <_R2.fq.gz> --sample <OUTPUT_SAMPLE_NAME> [--outdir <OUTDIR> --adapter <ADAPTER.fa> --cores <NCORES> --star <STAR_path> --genome_dir <STAR index path>]"
}

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --r1) # forward.fastq.gz, required
            R1="$2"
            shift
            shift
            ;;
        --r2) # reverse.fastq.gz, required
            R2="$2"
            shift
            shift
            ;;
        --sample) # output sample name, required
            SAMPLE="$2"
            shift
            shift
            ;;
        --outdir) # output directory, optional, default 'bams'
            OUTDIR="$2"
            shift
            shift
            ;;
        --adapter) # adapter fasta file, default to TruSeq3-PE-2.fa
            ADAPTER="$2"
            shift
            shift
            ;;
        --cores) # number of cores. default to 16.
            CORES="$2"
            shift
            shift
            ;;
        --star) # file path for STAR aligner.
            STAR="$2"
            shift
            shift
            ;;
        --genome_dir) # path to genome index directory.
            GENOME_DIR="$2"
            shift
            shift
            ;;
        *)
            # Unknown parameter
            echo "Unknown parameter: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$SAMPLE" ]; then
    echo "Missing required parameters!"
    show_usage
    exit 1
fi

########
# trim #
########
mkdir ${OUTDIR}

java -jar /home/yuanh/programs/source/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${CORES} \
    ${R1} ${R2} \
    ${OUTDIR}/${SAMPLE}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE}_R1_unpaired.fq.gz \
    ${OUTDIR}/${SAMPLE}_R2_paired.fq.gz ${OUTDIR}/${SAMPLE}_R2_unpaired.fq.gz \
    ILLUMINACLIP:${ADAPTER}:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

##########
# fastqc #
##########
fastqc -q -t ${CORES} -o ${OUTDIR} ${OUTDIR}/${SAMPLE}_R1_paired.fq.gz
fastqc -q -t ${CORES} -o ${OUTDIR} ${OUTDIR}/${SAMPLE}_R2_paired.fq.gz

#############
# alignment #
#############
# params based on ENCODE: https://github.com/ENCODE-DCC/rna-seq-pipeline/blob/dev/src/align.py

mkdir ${STAR_out}
${STAR} --genomeDir ${GENOME_DIR} \
    --readFilesIn ${OUTDIR}/${SAMPLE}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE}_R2_paired.fq.gz \
    --readFilesCommand zcat \
    --runThreadN ${CORES} \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverReadLmax 1 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterType BySJout \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbScore 1 \
    --limitBAMsortRAM ${ramGB}000000000 \
    --outStd BAM_SortedByCoordinate \
    --outFileNamePrefix ${STAR_out}/${SAMPLE}_ > ${OUTDIR}/${SAMPLE}.bam

##########
# filter #
##########
# filtering:
# -f 0x02 only keep reads mapped in proper pair. 
# -q 20 only keep fragments where the total q of two mates are > 20.

BAM=${OUTDIR}/${SAMPLE}.bam
samtools view -h ${BAM}|samtools view -bh -f 0x02 -q 20 -S - > ${BAM/.bam/.filter.bam}
samtools index ${BAM/.bam/.filter.bam}
samtools flagstat ${BAM/.bam/.filter.bam} > ${BAM/.bam/.filter.flagstat}

# remove intermediate files
rm ${BAM}
rm ${OUTDIR}/${SAMPLE}_*.gz