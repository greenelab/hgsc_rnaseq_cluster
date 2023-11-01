require("tximeta")
require("data.table")
require("dplyr")

# this script reads in all SALMON quant files
# creates a SummarizedExperiment object
# adds in the metadata
# summarizes transcript abundance to gene abundance
# this code is heavily inspired by
# https://github.com/AlexsLemonade/training-modules/blob/master/RNA-seq/02-gastric_cancer_tximeta.Rmd


#' Read in the Quant files
#'
#' @param quant_dir The directory containing all salmon quant files
#' the files should be within a directory titled by the sample id
#' @return a dataframe of the sample name and full file path
read_quant_files <- function(quant_dir){

    res_files = list.files(quant_dir,
                            pattern=".sf",
                            full.names = TRUE,
                            recursive=T)

    # split on the path, we assume it is the name of the folder the sf file is located
    sample_names = unlist(lapply(strsplit(res_files, "/"),
                            function(x) x[[length(x)-1]]))

    # now this name is too verbose, split on _ and take the first element
    # and the seventh element (the lane number)
    sample_date = unlist(lapply(strsplit(sample_names, "_"),
                            function(x) x[[2]]))
    sample_names = unlist(lapply(strsplit(sample_names, "_"),
                            function(x) paste("Sample", x[[1]], x[[7]], sep="_")))

    # now make the dataframe
    tximeta_df = data.frame(files = res_files,
                            names = sample_names,
                            date = sample_date)
    return(tximeta_df)
}


#' Transcript to Genes
#'
#' @param metadata_df The file with all needed metadata
#' to understand the counts
#' @param a dataframe of the sample name and full file path for tximeta
#' @return a Summarizedexperiment of gene abundance estimation
process_tx2gene <- function(metadata_df, quant_df){

    # first we merge together the column data
    metadata_df$lane_id = paste0("L00", metadata_df$platform_unit_id)
    metadata_df$names = paste("Sample", metadata_df$sample_id,
                                metadata_df$lane_id, sep="_")
    metadata_df = unique( metadata_df[,c("names", "experimental_strategy", "platform")])
    col_data = dplyr::inner_join(metadata_df, quant_df, by = "names")

    # now read in
    txi_data <- tximeta(col_data)

    # summarize
    gene_summarized <- summarizeToGene(txi_data)

    return(gene_summarized)
}


### read in arguements
args = commandArgs(trailingOnly=TRUE)
quant_dir = args[1]
metadata_file = args[2]
outdir = args[3]

metadata_df <- data.frame(fread(metadata_file))
quant_df <- read_quant_files(quant_dir)
gene_summarized <- process_tx2gene(metadata_df, quant_df)

txi_out_file = file.path(outdir, "salmon_gene_quant.RDS")
readr::write_rds(gene_summarized, file = txi_out_file)

# write it out in chunks for upload to github
for(idx in 1:9){
    file_name = paste0("salmon_gene_quant_", idx, ".RDS")
    txi_out_file = file.path(outdir, file_name)

    start = (idx-1)*35+1
    end = idx*35
    readr::write_rds(gene_summarized[,start:end], file = txi_out_file)

}
