
#' Get translation table for gene ids of interest
#'
#' @param gene_vec, vector of gene ids to translate
#' @param filter_type, type of id used in genevec
#' @param attributes, the gene ids you are interested in returning
#' @return gene_ids_df, data frame of ensembl
#' and hgnc ids filtered by gene_vec
get_gene_id_map <- function(gene_vec, filter_type, attributes, all_chr=F){
    require('biomaRt')

    # must include filter type, in attributes
    # and chromosome
    attributes = na.omit(unique(c(filter_type, "chromosome_name", attributes)))

    # translate the gene names
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="hsapiens_gene_ensembl" ,
                    host="https://useast.ensembl.org")

    mart <- useDataset("hsapiens_gene_ensembl",
                    ensembl)

    if(is.na(filter_type)){
        gene_ids_df <- getBM(attributes= attributes,
                             mart= mart)

    }else{
        gene_ids_df <- getBM(filters= filter_type,
                             attributes= attributes,
                             values=gene_vec,
                             mart= mart)
    }



    # remove things from non-typical chromosomes
    if(!all_chr){
        # get all the genes with a non-canonical chromosome
        non_reg_chr = grep("HG|HS|NOVEL|PATCH|GL", unique(gene_ids_df$chromosome_name), value=T, invert=F)
        gene_ids_df = gene_ids_df[which(!gene_ids_df$chromosome_name %in% non_reg_chr),]

    }else if(length(intersect(c("ensembl_gene_id", "hgnc_symbol"), attributes) == 2)){

        # only remove from non-typical chromosome
        # if the gene is on another chromosome

        # remove genes with no hgnc_symbol
        gene_ids_df = subset(gene_ids_df, hgnc_symbol != "")

        # get all the genes with a non-canonical chromosome
        non_reg_chr = grep("HG|HS|NOVEL|PATCH|GL", unique(gene_ids_df$chromosome_name), value=T, invert=F)
        hgnc_nonreg = gene_ids_df[which(gene_ids_df$chromosome_name %in% non_reg_chr),]
        hgnc_reg = gene_ids_df[which(!gene_ids_df$chromosome_name %in% non_reg_chr),]

        # now only keep noncanonical IF it only has nonreg
        hgnc_has_only_nonreg = setdiff(unique(hgnc_nonreg$hgnc_symbol), unique(hgnc_reg$hgnc_symbol))
        hgnc_nonreg = hgnc_nonreg[which(hgnc_nonreg$hgnc_symbol %in% hgnc_has_only_nonreg),]

        # for the nonreg keep the most updated one
        hgnc_nonreg = hgnc_nonreg[order(hgnc_nonreg$ensembl_gene_id), ]
        hgnc_nonreg_dup = which(!duplicated(hgnc_nonreg$hgnc_symbol))
        hgnc_nonreg = hgnc_nonreg[hgnc_nonreg_dup, ]

        gene_ids_df = rbind(hgnc_reg, hgnc_nonreg)
    }

    return(gene_ids_df)

}

#' Read and format nanostring data
#'
#' @param nano_expr_file, file path to nanostring data
#' @return melted and translated data frame
#' with expression, suid, hgnc_symbol, and ensembl_id
format_nanostring_data <- function(nano_expr_file){

    require(data.table)

    # read in file
    nano_df = fread(nano_expr_file)
    colnames(nano_df)[1] = "suid"

    # melt table for easier processing
    nano_df_melt = melt(nano_df, id.vars="suid")
    colnames(nano_df_melt)[2:3] = c("hgnc_symbol", "expr")
    nano_df_melt$hgnc_symbol = toupper(nano_df_melt$hgnc_symbol)

    # translate gene ids
    gene_ids_df = get_gene_id_map(unique(nano_df_melt$hgnc_symbol),
                                    filter_type="hgnc_symbol",
                                    attributes= c("ensembl_gene_id", "hgnc_symbol"))
    nano_df_melt_trans = merge(gene_ids_df, nano_df_melt, all.y=T)

    return(nano_df_melt_trans[,c("suid", "ensembl_gene_id", "hgnc_symbol", "expr")])
}


#' Read and format hta data
#'
#' @param hta_expr_file, file path to hta data
#' @param hta_trans_file, file path to sample ID translation table
#' @return melted and translated data frame
#' with expression, suid, hgnc_symbol, and ensembl_id
format_hta_data <- function(hta_expr_file, hta_trans_file){

    require(data.table)

    # read in expression
    hta_df = fread(hta_expr_file)
    colnames(hta_df)[1] = "ensembl_gene_id"

    #  read sample translation file
    hta_trans_df = fread(hta_trans_file)
    colnames(hta_trans_df)[1:2] = c("name", "suid")
    hta_trans_df = hta_trans_df[,1:2]

    # melt the table to easier processing
    hta_df_melt = melt(hta_df, id.vars="ensembl_gene_id")
    colnames(hta_df_melt)[2:3] = c("name", "expr")

    # translate gene ids
    # first we need to remove the _at in the ensembl_gene_ids
    hta_df_melt$ensembl_gene_id = gsub("_at", "", hta_df_melt$ensembl_gene_id)
    gene_ids_df = get_gene_id_map(unique(hta_df_melt$ensembl_gene_id),
                                    filter_type="ensembl_gene_id",
                                    attributes= c("hgnc_symbol", "ensembl_gene_id"))
    hta_df_melt_trans = merge(as.data.table(gene_ids_df), hta_df_melt, all.y=T, by="ensembl_gene_id")

    # translate sample ids
    hta_df_melt_trans_trans = merge(hta_trans_df, hta_df_melt_trans, all.y=T, by="name")

    return(hta_df_melt_trans_trans[,c("suid", "ensembl_gene_id", "hgnc_symbol", "expr")])
}


#' Read and format RNASeq data
#'
#' @param rnaseq, file path to hta data
#' @param rnaseq_trans_file, file path to sample ID translation table
#' @return melted and translated data frame
#' with expression, suid, hgnc_symbol, and ensembl_id
format_rnaseq_data <- function(rnaseq_expr_file, rnaseq_trans_file, extra_meta_df = NA, isNCOCS=F){ #### EDIT HERE

    library(readr)
    library(SummarizedExperiment)


    # read in the RDS
    gene_summarized = readr::read_rds(rnaseq_expr_file)
    gene_count = SummarizedExperiment::assay(gene_summarized, "counts")

    # read in translation file
    if(!is.na(rnaseq_trans_file)){
        rnaseq_samp_trans = fread(rnaseq_trans_file)
        if(isNCOCS){
            rnaseq_samp_trans = unique(rnaseq_samp_trans[,c("ID", "BLOCK_ID")])
            colnames(rnaseq_samp_trans) = c("ID", "suid")
            rnaseq_samp_trans = na.omit(rnaseq_samp_trans)
        }else{
            rnaseq_samp_trans = unique(rnaseq_samp_trans[,c("ID", "suid")])
        }

    }else{
        rnaseq_samp_trans = data.frame(ID=colnames(gene_count),
                                       suid=colnames(gene_count))
    }


    # filter lowly expressed genes (median across samples must be more than 0)
    filtered_genes = rownames(gene_count)[which(apply(gene_count, 1, median)==0)]
    gene_count_filt = subset(gene_count,
                            ! rownames(gene_count) %in% filtered_genes)

    # get the 85th quantile and library size normalize
    samp_85th = apply(gene_count_filt, 2, quantile, 0.85)
    samp_85th = samp_85th/max(samp_85th)
    gene_count_filt_norm = sweep(gene_count_filt, 2, samp_85th, FUN = '/')

    # metl for easier processing
    gene_count_filt_norm_melt = data.frame(gene_count_filt_norm)
    gene_count_filt_norm_melt$gene_id = row.names(gene_count_filt_norm_melt)
    gene_count_filt_norm_melt <- melt(gene_count_filt_norm_melt, id.vars="gene_id")
    colnames(gene_count_filt_norm_melt) = c("ensembl_gene_id", "ID", "expr")


    # translate gene ids
    # first we need to remove the _at in the ensembl_gene_ids
    gene_ids_df = get_gene_id_map(unique(gene_count_filt_norm_melt$ensembl_gene_id),
                                    filter_type="ensembl_gene_id",
                                    attributes= c("hgnc_symbol", "ensembl_gene_id"))
    gene_count_filt_norm_melt_trans = merge(as.data.table(gene_ids_df), gene_count_filt_norm_melt, all.y=T, by="ensembl_gene_id")

    # translate sample ids
    gene_count_filt_norm_melt_trans$ID = gsub("Sample_", "", gene_count_filt_norm_melt_trans$ID)

    if(isNCOCS){
        rnaseq_samp_trans$ID = as.character(rnaseq_samp_trans$ID)
        rnaseq_samp_trans$suid = as.character(rnaseq_samp_trans$suid)
    }
    final_dt = merge(rnaseq_samp_trans, gene_count_filt_norm_melt_trans, all.y=T, by="ID")

    if(is.na(rnaseq_trans_file)){
        final_dt$suid = final_dt$ID
    }

    colnames_return = c("suid", "ID", "ensembl_gene_id", "hgnc_symbol", "expr")
    if(!is.na(extra_meta_df)){
        final_dt = merge(final_dt, extra_meta_df, all.x=T, by="ID")
        colnames_return = c(colnames_return, setdiff(colnames(extra_meta_df), "ID"))
    }

    return(final_dt[,..colnames_return])
}
