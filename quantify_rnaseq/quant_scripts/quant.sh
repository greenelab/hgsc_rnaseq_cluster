#!/bin/bash

#### This file runs fastp which also performs trimming

IN_DIR=$1
INDEX_DIR=$2
OUT_DIR=$3
SAMPLE_ID=$4

MATE_1=${IN_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz
MATE_2=${IN_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz

# setup environment
module load anaconda3/3.16.0
source activate aaces_quant

salmon quant -i ${INDEX_DIR} \
    -l A \
    -1 ${MATE_1} \
    -2 ${MATE_2} \
    -o ${OUT_DIR} \
    --validateMappings --rangeFactorizationBins 4 \
    --gcBias --seqBias \
    --threads 4