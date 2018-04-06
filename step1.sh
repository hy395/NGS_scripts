#!/bin/bash

## initial fastqc
## this is updated on 4.18.2017
## Arguments:
## 1. input directory
## 2. output directory

INDIR=$1
OUTDIR=$2

mkdir ${OUTDIR}
for ENTRY in ${INDIR}/*.gz
do
 fastqc -q -t 8 -o ${OUTDIR} ${ENTRY}
 echo ${ENTRY}
done

