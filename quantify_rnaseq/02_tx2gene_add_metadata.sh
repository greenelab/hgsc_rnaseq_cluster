#!/bin/bash

## this script converts transcript to gene abundance
## and adds metadata information
## 2 files are written
## 1) gene abundance (RDs)
## 2) transcript abundance (RDS)

ROOT_DIR=$1
RUN_FULL=$2

# ROOT_DIR="/YOUR_PATH/hgsc_characterization/"
# RUN_FULL=False means that only preliminary data will be run


RNASEQ_DIRNAME="rna_seq_full"
MANIFEST_FILENAME="19291R_MANIFEST"

if [ "$RUN_FULL" = "False" ]; then
    echo "only running the preliminary data"
    RNASEQ_DIRNAME="rna_seq"
    MANIFEST_FILENAME="18341R_MANIFEST"
else
    echo "Running full data"
fi

##################
## set up directories
##################

QUANT_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/salmon_quant/
METADATA_FILE=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/${MANIFEST_FILENAME}.csv
OUT_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/salmon_quant_processed/

if [ ! -d ${OUT_DIR} ]; then
    mkdir -p ${OUT_DIR}
fi


##################
## process and write out files
##################
module load R/4.0.3
process_script=${ROOT_DIR}/quantify_rnaseq/tx2gene_add_metadata_scripts/tximeta_process.R
Rscript ${process_script} ${QUANT_DIR} ${METADATA_FILE} ${OUT_DIR}