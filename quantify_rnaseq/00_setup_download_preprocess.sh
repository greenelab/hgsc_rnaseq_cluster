#!/bin/bash

## this script sets up everything for salmon quantification
## 1) downloads and validates the raw fastq files
## 2) runs fastp for Qc and trimming
## 3) downloads HG38 reference files
## 3) builds salmon index
## it also assumes you will be running this on an lsf cluster

ROOT_DIR=$1
LOG_DIR=$2
RUN_FULL=$3

# cd /YOUR_PATH/hgsc_characterization/data/rna_seq_full/fastq_files
# ROOT_DIR="/YOUR_PATH/hgsc_characterization/"
# LOG_DIR="/YOUR_PATH/logs/"
# RUN_FULL=False means that only preliminary data will be run

RNASEQ_DIRNAME="rna_seq_full"
WGET_FILENAME="wget_cmd_rna_full3"
MANIFEST_FILENAME="19291R_MANIFEST"
if [ "$RUN_FULL" = "False" ]; then
    echo "only running the preliminary data"
    RNASEQ_DIRNAME="rna_seq"
    WGET_FILE="wget_cmd_rna_trial"
    MANIFEST_FILENAME="18341R_MANIFEST"
else
    echo "Running full data"
fi

##################
## download all rna-seq files
################## 

sh ${ROOT_DIR}/quantify_rnaseq/setup_download_preprocess_scripts/${WGET_FILENAME}.sh
md5sum ${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/*.fastq.gz > ${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/md5sums.txt
module load R/4.0.3
Rscript ${ROOT_DIR}/quantify_rnaseq/setup_download_preprocess_scripts/check_md5sum.R \
        ${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/${MANIFEST_FILENAME}.csv \
        ${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/md5sums.txt



##################
## set up directories
##################

IN_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_files/
TRIM_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastq_trimmed/
QC_DIR=${ROOT_DIR}/data/${RNASEQ_DIRNAME}/fastp_out/

##################
## run fastp on all raw files
##################

fastp_script=${ROOT_DIR}/quantify_rnaseq/setup_download_preprocess_scripts/pre_quant_qc.sh
for fastq_file in ${IN_DIR}/*R1*.fastq.gz
do
    sample_id=$(basename ${fastq_file})
    sample_id=${sample_id:0:${#sample_id}-16}

    # only run it if we haven't done it before
    JSON_OUT=${QC_DIR}/${SAMPLE_ID}_fastp.json
    if [ ! -s ${JSON_OUT} ]; then
        echo running: $sample_id

        bsub -R "select[mem>10] rusage[mem=10]" -J ${sample_id} \
                -o ${LOG_DIR}/stdout_%J.out -e ${LOG_DIR}/stderr_%J.err \
                ${fastp_script} \
                ${IN_DIR} ${TRIM_DIR} ${QC_DIR} ${sample_id}
    fi

done

# make a summary of all fastp files
module load anaconda3/3.16.0
source activate aaces_quant

cd ${QC_DIR}
multiqc .
cd -

##################
# download references to make index
# download scripts heavily inspired by
# https://github.com/AlexsLemonade/training-txome-prep
##################

## download reference data if not available
# first fasta
cdna_dir=${ROOT_DIR}/reference_data/cdna
mkdir -p ${cdna_dir}

if [ ! -s ${cdna_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz ]; then
    wget --directory-prefix=${cdna_dir} ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
fi

# now get gtfs
gtf_dir=${ROOT_DIR}/reference_data/gtf
mkdir ${gtf_dir}
if [ ! -s ${gtf_dir}/Homo_sapiens.GRCh38.95.gtf ]; then
    wget --directory-prefix=${gtf_dir} ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz
    gunzip ${gtf_dir}/*
fi

## now if both succeeded build the index

index_dir=${ROOT_DIR}/reference_data/index/Homo_sapiens
mkdir -p ${index_dir}

if [ ! -s ${index_dir}/long_index/ref_k31_fixed.fa ] && [ -s ${cdna_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz ] ; then
    # use long kmers
    salmon index \
        -t ${cdna_dir}/Homo_sapiens.GRCh38.cdna.all.fa.gz \
        -i ${index_dir}/long_index \
        -k 31
fi

## now run the tx2gene
tx2gene_dir=${ROOT_DIR}/reference_data/tx2gene
mkdir -p ${tx2gene_dir}

Rscript ${ROOT_DIR}/quantify_rnaseq/setup_download_preprocess_scripts/get_tx2gene.R \
  --gtf_file ${gtf_dir}/Homo_sapiens.GRCh38.95.gtf \
  --output_file ${tx2gene_dir}/Homo_sapiens.GRCh38.95_tx2gene.tsv

