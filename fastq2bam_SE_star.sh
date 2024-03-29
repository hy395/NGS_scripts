#!/bin/bash

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

#############
# arguments #
#############

# param defaults
OUTDIR='bams'
ADAPTER='/home/yuanh/programs/source/Trimmomatic-0.39/adapters/TruSeq3-SE.fa' # default for Illumina RNA-seq
CORES=16
GENOME_DIR='/scratch4/yuanh/star_index/GRCh38_Ensembl99_sparseD3_sjdbOverhang99'
STAR='/home/yuanh/programs/source/STAR-2.7.3a/bin/Linux_x86_64_static/STAR'
STAR_out='STAR_out'
ramGB=20

function show_usage() {
    echo "Usage: ./fastq2bam_SE_star.sh --fastq <fq.gz> --sample <OUTPUT_SAMPLE_NAME> [--outdir <OUTDIR> --adapter <ADAPTER.fa> --cores <NCORES> --star <STAR_path> --genome_dir <STAR index path>]"
}


while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fastq) # _.fastq.gz, required
            FQ="$2"
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
        --adapter) # adapter fasta file, default to TruSeq3-SE.fa
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
if [ -z "$FQ" ] || [ -z "$SAMPLE" ]; then
    echo "Missing required parameters!"
    show_usage
    exit 1
fi

# trim
java -jar /home/yuanh/programs/source/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads ${CORES} \
    ${FQ} \
    ${OUTDIR}/${SAMPLE}_trim.fastq.gz \
    ILLUMINACLIP:${ADAPTER}:2:30:7 SLIDINGWINDOW:3:15 MINLEN:36

##########
# fastqc #
##########
fastqc -q -t ${CORES} -o ${OUTDIR} ${OUTDIR}/${SAMPLE}_trim.fastq.gz

#############
# alignment #
#############
# params based on ENCODE: https://github.com/ENCODE-DCC/rna-seq-pipeline/blob/dev/src/align.py
${STAR} --genomeLoad NoSharedMemory --genomeDir ${GENOME_DIR} \
    --readFilesIn ${OUTDIR}/${SAMPLE}_trim.fastq.gz \
    --readFilesCommand zcat \
    --runThreadN ${CORES} \
    --outFileNamePrefix ${STAR_out}/${SAMPLE}_ \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterType BySJout \
    --sjdbScore 1 \
    --limitBAMsortRAM ${ramGB}000000000
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate \
    --outSAMstrandField intronMotif > ${OUTDIR}/${SAMPLE}.bam

#############
# filtering #
#############
# filtering: -F 0x204 removes 1) reads that don't map. 2) QC failed reads.
BAM=${OUTDIR}/${SAMPLE}.bam
samtools view -h ${BAM}|samtools view -bh -F 0x204 -q 10 -S - > ${BAM/.bam/.filter.bam}
samtools index ${BAM/.bam/.filter.bam}
samtools flagstat ${BAM/.bam/.filter.bam} > ${BAM/.bam/.filter.flagstat}

rm ${BAM}
rm ${OUTDIR}/${SAMPLE}_trim.fastq.gz
