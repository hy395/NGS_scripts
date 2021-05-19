#!/bin/bash

# step1. concatenate all fastq files
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_1_*_R1_* | gzip > raw/ATAC_4.10_A_1_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_1_*_R2_* | gzip > raw/ATAC_4.10_A_1_R2.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_2_*_R1_* | gzip > raw/ATAC_4.10_A_2_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_2_*_R2_* | gzip > raw/ATAC_4.10_A_2_R2.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_3_*_R1_* | gzip > raw/ATAC_4.10_A_3_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ Adherent\ ATACseq/AML_4_10Adherent_3_*_R2_* | gzip > raw/ATAC_4.10_A_3_R2.fastq.gz &
#
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_1_*_R1_* | gzip > raw/ATAC_4.10_S_1_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_1_*_R2_* | gzip > raw/ATAC_4.10_S_1_R2.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_2_*_R1_* | gzip > raw/ATAC_4.10_S_2_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_2_*_R2_* | gzip > raw/ATAC_4.10_S_2_R2.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_3_*_R1_* | gzip > raw/ATAC_4.10_S_3_R1.fastq.gz &
#zcat Papapetrou/ATACseq/AML-4.10\ suspension\ ATACseq/AML_4_10Suspension_3_*_R2_* | gzip > raw/ATAC_4.10_S_3_R2.fastq.gz &
#
#
#zcat Papapetrou/RNAseq/AML-4.10\ adherent\ RNAseq/AML_4_10Adherent_1_*_R1_* | gzip > raw/RNA_4.10_A_1_R1.fastq.gz &
#zcat Papapetrou/RNAseq/AML-4.10\ adherent\ RNAseq/AML_4_10Adherent_2_*_R1_* | gzip > raw/RNA_4.10_A_2_R1.fastq.gz &
#zcat Papapetrou/RNAseq/AML-4.10\ adherent\ RNAseq/AML_4_10Adherent_3_*_R1_* | gzip > raw/RNA_4.10_A_3_R1.fastq.gz &
#zcat Papapetrou/RNAseq/AML-4.10\ suspension\ RNAseq/AML_4_10Suspension_1_*_R1_* | gzip > raw/RNA_4.10_S_1_R1.fastq.gz &
#zcat Papapetrou/RNAseq/AML-4.10\ suspension\ RNAseq/AML_4_10Suspension_2_*_R1_* | gzip > raw/RNA_4.10_S_2_R1.fastq.gz &
#zcat Papapetrou/RNAseq/AML-4.10\ suspension\ RNAseq/AML_4_10Suspension_3_*_R1_* | gzip > raw/RNA_4.10_S_3_R1.fastq.gz &

# step2. run fastqc
#~/programs/myscripts/step1.sh raw raw_fastqc

############
# ATAC-seq #
############
# step3: trim, align, filter and rmdup for ATAC PE data.
#mkdir ATAC_clean_bam
#for i in ATAC_4.10_A_1 ATAC_4.10_A_2 ATAC_4.10_A_3 ATAC_4.10_S_1 ATAC_4.10_S_2 ATAC_4.10_S_3; do
# /home/hy395/programs/myscripts/fastq2bam_PE_bowtie.sh /data/leslie/genome_index/hg19.bowtie2/hg19 raw/${i}_R1.fastq.gz raw/${i}_R2.fastq.gz /home/hy395/programs/source/Trimmomatic-0.36/adapters/NexteraPE-PE.fa ATAC_clean_bam/${i} 20 > ATAC_clean_bam/${i}.log 2>&1
# echo ${i}
#done

# step4: peak calling and IDR
# 1) convert to bed, shift
#mkdir ATAC_shift_beds
#for i in ATAC_4.10_A_1 ATAC_4.10_A_2 ATAC_4.10_A_3 ATAC_4.10_S_1 ATAC_4.10_S_2 ATAC_4.10_S_3; do
# samtools sort -n -T ${i}_nsorted ATAC_clean_bam/${i}.filter.rmdup.bam |bedtools bamtobed -i - | /home/hy395/programs/myscripts/adjustBedTn5.sh - > ATAC_shift_beds/${i}_SE.bed
# echo ${i}
#done

# 2) peak calling
#mkdir macs2out
#cat ATAC_shift_beds/ATAC_4.10_A_*_SE.bed > ATAC_shift_beds/ATAC_4.10_A_pool.bed
#cat ATAC_shift_beds/ATAC_4.10_S_*_SE.bed > ATAC_shift_beds/ATAC_4.10_S_pool.bed
#macs2 callpeak -t ATAC_shift_beds/ATAC_4.10_A_pool.bed -f BED -n A_pool --outdir macs2out --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
#macs2 callpeak -t ATAC_shift_beds/ATAC_4.10_S_pool.bed -f BED -n S_pool --outdir macs2out --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
#for i in ATAC_4.10_A_1 ATAC_4.10_A_2 ATAC_4.10_A_3 ATAC_4.10_S_1 ATAC_4.10_S_2 ATAC_4.10_S_3; do
# macs2 callpeak -t ATAC_shift_beds/${i}_SE.bed -f BED -n ${i} --outdir macs2out --keep-dup all -g hs -p 1e-2 --nomodel --shift -37 --extsize 73
#done

# 3) IDR
#idr --samples macs2out/ATAC_4.10_A_1_peaks.narrowPeak macs2out/ATAC_4.10_A_2_peaks.narrowPeak --peak-list macs2out/A_pool_peaks.narrowPeak --output-file idr_A_1_2.txt
#idr --samples macs2out/ATAC_4.10_A_1_peaks.narrowPeak macs2out/ATAC_4.10_A_3_peaks.narrowPeak --peak-list macs2out/A_pool_peaks.narrowPeak --output-file idr_A_1_3.txt
#idr --samples macs2out/ATAC_4.10_A_2_peaks.narrowPeak macs2out/ATAC_4.10_A_3_peaks.narrowPeak --peak-list macs2out/A_pool_peaks.narrowPeak --output-file idr_A_2_3.txt
#idr --samples macs2out/ATAC_4.10_S_1_peaks.narrowPeak macs2out/ATAC_4.10_S_2_peaks.narrowPeak --peak-list macs2out/S_pool_peaks.narrowPeak --output-file idr_S_1_2.txt
#idr --samples macs2out/ATAC_4.10_S_1_peaks.narrowPeak macs2out/ATAC_4.10_S_3_peaks.narrowPeak --peak-list macs2out/S_pool_peaks.narrowPeak --output-file idr_S_1_3.txt
#idr --samples macs2out/ATAC_4.10_S_2_peaks.narrowPeak macs2out/ATAC_4.10_S_3_peaks.narrowPeak --peak-list macs2out/S_pool_peaks.narrowPeak --output-file idr_S_2_3.txt

#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_A_1_2.txt hg19 0.05
#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_A_1_3.txt hg19 0.05
#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_A_2_3.txt hg19 0.05
#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_S_1_2.txt hg19 0.05
#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_S_1_3.txt hg19 0.05
#Rscript ~/programs/myscripts/filter.R ATAC_idr/idr_S_2_3.txt hg19 0.05

# additional. generate shifted and centered bw and rds files
for i in ATAC_4.10_A_1 ATAC_4.10_A_2 ATAC_4.10_A_3 ATAC_4.10_S_1 ATAC_4.10_S_2 ATAC_4.10_S_3; do
 Rscript ~/programs/myscripts/makeBW_atac.R ATAC_clean_bam/${i}.filter.rmdup.bam 73
done

###########
# RNA-seq #
###########
# step3: trim, align for RNA-seq SE data. there is no adapter detected in fastq, so adapter.fa is empty file.
#mkdir star_alignment
#mkdir RNA_clean_bam
#for i in RNA_4.10_A_1_R1 RNA_4.10_A_2_R1 RNA_4.10_A_3_R1 RNA_4.10_S_1_R1 RNA_4.10_S_2_R1 RNA_4.10_S_3_R1; do
# /home/hy395/programs/myscripts/fastq2bam_SE_star.sh raw/${i}.fastq.gz ${i} RNA_clean_bam adapter.fa /data/leslie/genome_index/hg19.star star_alignment 1 > RNA_clean_bam/${i}.log 2>&1
# echo ${i}
#done

# step4: count
#Rscript count.R

#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_A_1_R1.filter.bam 0
#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_A_2_R1.filter.bam 0
#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_A_3_R1.filter.bam 0
#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_S_1_R1.filter.bam 0
#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_S_2_R1.filter.bam 0
#Rscript ~/programs/myscripts/makeBW.R RNA_clean_bam/RNA_4.10_S_3_R1.filter.bam 0


