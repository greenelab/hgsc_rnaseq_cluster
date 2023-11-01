# This file was copied into this repo on March 11, 2021
# from https://github.com/AlexsLemonade/training-txome-prep

# J. Taroni for CCDL 2019
# Purpose: generating a tx2gene TSV file for use with tximport::tximport from
# an Ensembl GTF file
#
# USAGE: Rscript scripts/get_tx2gene.R \
#         --gtf_file <PATH_TO_ENSEMBL_GTF_FILE> \
#         --output_file <PATH_TO_OUTPUT_TSV_FILE>
#
#
# Note this will generate a large .sqlite file, the size of which depends on
# the species

#### install optparse if necessary ---------------------------------------------

if (!("optparse" %in% installed.packages())) {
    install.packages("optparse", suppressUpdates = TRUE)
}

#### option parsing ------------------------------------------------------------

option_list <- list(
  optparse::make_option("--gtf_file", type = "character"),
  optparse::make_option("--output_file", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
gtf_file <- opt$gtf_file
output_file <- opt$output_file

#### main ----------------------------------------------------------------------

db_file <- sub("\\.gtf", ".sqlite", gtf_file)
if (!file.exists(db_file)) {
  ensembldb::ensDbFromGtf(gtf_file, outfile = db_file)
}
edb <- ensembldb::EnsDb(db_file)
tx <- ensembldb::transcriptsBy(edb, by = "gene")
tx.df <- as.data.frame(tx@unlistData)
tx2gene <- tx.df[, c("tx_name", "gene_id")]
readr::write_tsv(tx2gene, path = output_file)
