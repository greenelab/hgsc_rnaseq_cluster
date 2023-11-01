#!/bin/bash

## this script quantifies transcript expression using salmon
## it assumes 00_setup_download_preprocess.sh has already been run successfully
## is also assumes you will be running this on an lsf cluster
ROOT_DIR=$1
LOG_DIR=$2
RUN_FULL=$3

# ROOT_DIR="/YOUR_PATH/hgsc_characterization/"
# LOG_DIR="/YOUR_PATH/logs/"
# RUN_FULL=False means that only preliminary data will be run


RNASEQ_DIRNAME="rna_seq_full"
if [ "$RUN_FULL" = "False" ]; then
    echo "only running the preliminary data"
    RNASEQ_DIRNAME="rna_seq"
else
    echo "Running full data"
fi

##################
## set up directories
##################

IN_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_trimmed/
INDEX_DIR=${ROOT_DIR}/reference_data/index/Homo_sapiens/long_index/
OUT_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/salmon_quant/


##################
## run salmon on all trimmed fastq files
##################

quant_script=${ROOT_DIR}/quantify_rnaseq/quant_scripts/quant.sh
for fastq_file in ${TRIM_DIR}/*R1*.fastq.gz
do
    sample_id=$(basename ${fastq_file})
    sample_id=${sample_id:0:${#sample_id}-20}

    echo running: ${sample_id}

    CURR_OUT_DIR=${OUT_DIR}/${sample_id}
    if [ ! -d ${CURR_OUT_DIR} ]; then
        mkdir -p ${CURR_OUT_DIR}
    fi


    bsub -R "select[mem>10] rusage[mem=10]" -J ${sample_id} \
            -o ${LOG_DIR}/stdout_%J.out -e ${LOG_DIR}/stderr_%J.err \
            ${quant_script} \
            ${IN_DIR} ${INDEX_DIR} ${CURR_OUT_DIR} ${sample_id}

done
