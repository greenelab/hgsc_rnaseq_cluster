#!/bin/bash

#### This file runs fastp which also performs trimming

IN_DIR=$1
TRIM_DIR=$2
QC_DIR=$3
SAMPLE_ID=$4

MATE_1=${IN_DIR}/${SAMPLE_ID}_R1_001.fastq.gz
MATE_2=${IN_DIR}/${SAMPLE_ID}_R2_001.fastq.gz

TR_MATE_1=${TRIM_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz
TR_MATE_2=${TRIM_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz

JSON_OUT=${QC_DIR}/${SAMPLE_ID}_fastp.json
HTML_OUT=${QC_DIR}/${SAMPLE_ID}_fastp.html

echo ${MATE_1}
echo ${MATE_2}
echo ${TR_MATE_1}
echo ${TR_MATE_2}
echo ${JSON_OUT}
echo ${HTML_OUT}

# setup environment
module load anaconda3/3.16.0
source activate aaces_quant

# Run the adapter and quality trimming step -- also produces QC report
fastp -i ${MATE_1} \
    -I ${MATE_2} \
    -o ${TR_MATE_1} \
    -O ${TR_MATE_2} \
    --qualified_quality_phred 15 \
    --length_required 20 \
    --report_title ${SAMPLE_ID} \
    --json ${JSON_OUT} \
    --html ${HTML_OUT}