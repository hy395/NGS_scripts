# pre-generated index: https://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/
# wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/
# only compatible with STAR-2.7.1a

# or to generate:
STARDIR='/home/yuanh/programs/source/STAR-2.7.10b/bin/Linux_x86_64/STAR'
FASTA='/home/yuanh/programs/genomes/hg38/hg38.fa'
GTF='/home/yuanh/programs/genomes/hg38/gencode/v41/gencode.v41.basic.annotation.gtf'

${STARDIR} --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir ./ \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF}
