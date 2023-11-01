compare_centroids
================
Natalie Davidson
8/31/2023

# Running consensusOV

This notebook compares our clustering results with the output of
clusterOV.

## Write all methods out

``` r
# read microarray data
read_format_MA_expr <- function(in_df, metadata_table){
    
     rnaseq_expr_df = in_df

    # format it so we can add metadata
    gene_ids = row.names(rnaseq_expr_df)
    sample_ids = colnames(rnaseq_expr_df)
    rnaseq_expr_df = data.frame(t(rnaseq_expr_df))
    colnames(rnaseq_expr_df) = gene_ids
    
    # now add column names back so we can have annotations
    rnaseq_expr_df$ID = sample_ids
    
    full_df = merge(metadata_table, rnaseq_expr_df, by = "ID")
    
    return(list(full_df, gene_ids))
}

# read rnaseq
read_format_expr <- function(in_file, metadata_table){
    
     rnaseq_expr_df = data.frame(fread(in_file))

    # format it so we can add metadata
    gene_ids = rnaseq_expr_df[,1]
    rnaseq_expr_df = rnaseq_expr_df[,-1]
    sample_ids = colnames(rnaseq_expr_df)
    rnaseq_expr_df = data.frame(t(rnaseq_expr_df))
    colnames(rnaseq_expr_df) = gene_ids
    
    # now add column names back so we can have annotations
    rnaseq_expr_df$ID = gsub("Sample_", "", row.names(rnaseq_expr_df))
    
    full_df = merge(metadata_table, rnaseq_expr_df, by = "ID")
    
    return(list(full_df, gene_ids))
}

# filter scale RNASeq
filter_expr <- function(in_df, gene_df, metadata_table){
        
    # get samples that were used
    # get genes of interest
    gene_count_df = in_df
    gene_count_df = in_df[,gene_df$hgnc_symbol]
    gene_count_df = log10(gene_count_df+1)

    
    
    # add metadata back
    gene_count_df$ID = in_df$ID
    gene_count_df = merge(gene_count_df, metadata_table)
    
    return(gene_count_df)
}

# learn+match cluster labels between consensusOV and Way
get_cluster_labels <- function(gene_map, in_df, metadata_df, min_consensus_cutoff, do_log=FALSE){
    
    gene_map_curr = subset(gene_map, hgnc_symbol %in% gene_ids)
    gene_map_curr = unique(gene_map_curr[,c("hgnc_symbol", "entrezgene_id")])
    
    in_df = merge(in_df, gene_map_curr, by="hgnc_symbol" )
    
    # make numeric
    gene_ids_W = in_df$entrezgene_id
    in_df = subset(in_df, select = -c(hgnc_symbol, entrezgene_id))
    in_df = apply(as.matrix(in_df), 2, as.numeric)
    
    if(do_log){
        in_df = log10(in_df+1)
    }
    
    # do prediction
    consensus_res <- get.subtypes(in_df, gene_ids_W, method = "consensus")
    unstable_pred = names(which(apply(consensus_res$rf.probs, 1, max) < min_consensus_cutoff))
    
    # merge with old cluster labels
    res_df = data.frame(ID=colnames(in_df),
                             verhaak_clust=consensus_res$consensusOV.subtypes)
    res_df = merge(res_df, metadata_df)
    res_stable = subset(res_df, ! ID %in% unstable_pred)

    return(res_stable)
}


mapping_tcga_way <- function(x){
    if(x == "MES_consensus"){
        return(1)
    }else if(x == "PRO_consensus"){
        return(2)
    }else if(x == "IMR_consensus"){
        return(3)
    }else{
        return(4)
    }
}
```

## First read in all metadata for Schildkraut

``` r
# Black metadata
samp_meta_file = file.path(proj_dir, 
                    "/reference_data/main_AA_metadata_table.tsv")
metadata_table = data.frame(fread(samp_meta_file))

# format NAs for better plotting
metadata_table$failed_seq[is.na(metadata_table$failed_seq)] = "passed"
metadata_table$resequenced[is.na(metadata_table$resequenced)] = "passed"
metadata_table$low_qual[is.na(metadata_table$low_qual)] = "passed"
metadata_table$REMOVE_LOW_EXPRESSION[is.na(metadata_table$REMOVE_LOW_EXPRESSION)] = "passed"
metadata_table$REMOVE_LOW_EXPRESSION[metadata_table$REMOVE_LOW_EXPRESSION==TRUE] = "failed"

metadata_table_AA = metadata_table

# white metadata
samp_meta_file = file.path(proj_dir, 
                    "/reference_data/main_white_metadata_table.tsv")
metadata_table = data.frame(fread(samp_meta_file))

# format NAs for better plotting
metadata_table$failed_seq[is.na(metadata_table$failed_seq)] = "passed"
metadata_table$resequenced[is.na(metadata_table$resequenced)] = "passed"
metadata_table$low_qual[is.na(metadata_table$low_qual)] = "passed"
metadata_table$REMOVE_LOW_EXPRESSION[is.na(metadata_table$REMOVE_LOW_EXPRESSION)] = "passed"
metadata_table$REMOVE_LOW_EXPRESSION[metadata_table$REMOVE_LOW_EXPRESSION==TRUE] = "failed"

metadata_table_W = metadata_table


# genes of interest

common_genes_file = file.path(proj_dir, 
                                    "/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/CommonGenes_genelist.csv")
common_genes = data.frame(fread(common_genes_file))
colnames(common_genes) = "hgnc_symbol"


MAD_genes_file = file.path(proj_dir, 
                                    "/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/GlobalMAD_genelist.csv")
MAD_genes = data.frame(fread(MAD_genes_file))
colnames(MAD_genes) = "hgnc_symbol"

# cluster ID labels

new_file = file.path(proj_dir, 
                    "/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv")

clust_df = data.frame(fread(new_file))
colnames(clust_df)[1] = "ID"
```

## Read in expression for Schildkraut

``` r
# make formatted
in_file = file.path(proj_dir, 
                    "/data/rna_seq_pilot_and_new/salmon_normalized_filtered_for_way_pipeline.tsv")
res = read_format_expr(in_file, metadata_table_AA)
full_expr_B = res[[1]]

MAD_expr_B = filter_expr(full_expr_B, MAD_genes, metadata_table_AA)
MAD_expr_B = filter(MAD_expr_B, ran_in_way_pipeline == TRUE)
MAD_expr_B = subset(MAD_expr_B, select = -c(REMOVE_WHITE))
metadata_table_AA = subset(metadata_table_AA, select = -c(REMOVE_WHITE, REMOVE_NEOADJ))


# make formatted
in_file = file.path(proj_dir, 
                    "/data/rna_seq_whites/salmon_normalized_filtered_for_way_pipeline_whites.tsv")
res = read_format_expr(in_file, metadata_table_W)
full_expr_W = res[[1]]
MAD_expr_W = filter_expr(full_expr_W, MAD_genes, metadata_table_W)
MAD_expr_W = filter(MAD_expr_W, ran_in_way_pipeline == TRUE)
MAD_expr_W = subset(MAD_expr_W, select = -c(REMOVE_BLACK))
metadata_table_W = subset(metadata_table_W, select = -c(REMOVE_BLACK))

outlier_ids_W = MAD_expr_W$ID
```

# run clustOV

## run SchildkrautW

``` r
# format table
curr_IDs = MAD_expr_W$ID
gene_ids = colnames(MAD_expr_W)[2:ncol(MAD_expr_W)]
MAD_expr_W_matr = as.data.frame(t(MAD_expr_W[, 2:ncol(MAD_expr_W)]))
colnames(MAD_expr_W_matr) = curr_IDs

MAD_expr_W_matr$hgnc_symbol = gene_ids

# run clustering
schW_res = get_cluster_labels(gene_map, MAD_expr_W_matr, metadata_table_W, min_consensus_cutoff, do_log = TRUE)


schW_res$mapped = unlist(lapply(schW_res$verhaak_clust, mapping_tcga_way))
schW_consensus_way = schW_res$ID[which(schW_res$mapped == schW_res$ClusterK4_kmeans)]

table(schW_res$mapped, schW_res$ClusterK4_kmeans)
```

    ##    
    ##      1  2  3  4
    ##   1 50  0 12 10
    ##   2 12 29  1 11
    ##   3 39  3 34  9
    ##   4  5  8 37 49

``` r
# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(schW_res$mapped), 
                                reference = as.factor(schW_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 50  0 12 10
    ##          2 12 29  1 11
    ##          3 39  3 34  9
    ##          4  5  8 37 49
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.5243         
    ##                  95% CI : (0.467, 0.5811)
    ##     No Information Rate : 0.343          
    ##     P-Value [Acc > NIR] : 4.941e-11      
    ##                                          
    ##                   Kappa : 0.3581         
    ##                                          
    ##  Mcnemar's Test P-Value : 2.378e-08      
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.4717  0.72500   0.4048   0.6203
    ## Specificity            0.8916  0.91078   0.7733   0.7826
    ## Pos Pred Value         0.6944  0.54717   0.4000   0.4949
    ## Neg Pred Value         0.7637  0.95703   0.7768   0.8571
    ## Prevalence             0.3430  0.12945   0.2718   0.2557
    ## Detection Rate         0.1618  0.09385   0.1100   0.1586
    ## Detection Prevalence   0.2330  0.17152   0.2751   0.3204
    ## Balanced Accuracy      0.6817  0.81789   0.5890   0.7014

``` r
schW_ba = conf_res$byClass[,11]
```

## run SchildkrautB

``` r
# format table
curr_IDs = MAD_expr_B$ID
gene_ids = colnames(MAD_expr_B)[2:ncol(MAD_expr_B)]
MAD_expr_B_matr = as.data.frame(t(MAD_expr_B[, 2:ncol(MAD_expr_B)]))
colnames(MAD_expr_B_matr) = curr_IDs

MAD_expr_B_matr$hgnc_symbol = gene_ids

schB_res = get_cluster_labels(gene_map, MAD_expr_B_matr, metadata_table_AA, min_consensus_cutoff, do_log = TRUE)

table(schB_res$verhaak_clust, schB_res$ClusterK4_kmeans)
```

    ##                
    ##                  1  2  3  4
    ##   IMR_consensus 18  1 55  1
    ##   DIF_consensus 13 23 33 10
    ##   PRO_consensus  4 36  1  4
    ##   MES_consensus 44  3 12  4

``` r
schB_res$mapped = unlist(lapply(schB_res$verhaak_clust, mapping_tcga_way))
schB_consensus_way = schB_res$ID[which(schB_res$mapped == schB_res$ClusterK4_kmeans)]


# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(schB_res$mapped), 
                                reference = as.factor(schB_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 44  3 12  4
    ##          2  4 36  1  4
    ##          3 18  1 55  1
    ##          4 13 23 33 10
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.5534         
    ##                  95% CI : (0.491, 0.6146)
    ##     No Information Rate : 0.3855         
    ##     P-Value [Acc > NIR] : 2.761e-08      
    ##                                          
    ##                   Kappa : 0.4077         
    ##                                          
    ##  Mcnemar's Test P-Value : 5.666e-09      
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.5570   0.5714   0.5446  0.52632
    ## Specificity            0.8962   0.9548   0.8758  0.71605
    ## Pos Pred Value         0.6984   0.8000   0.7333  0.12658
    ## Neg Pred Value         0.8241   0.8756   0.7540  0.95082
    ## Prevalence             0.3015   0.2405   0.3855  0.07252
    ## Detection Rate         0.1679   0.1374   0.2099  0.03817
    ## Detection Prevalence   0.2405   0.1718   0.2863  0.30153
    ## Balanced Accuracy      0.7266   0.7631   0.7102  0.62118

``` r
schB_ba = conf_res$byClass[,11]
```

## run TCGA

``` r
tcga_data = data(TCGA_eset)
ExpressionData <- get(tcga_data)
tcga_dta <- exprs(ExpressionData)

tcga_metadata_table = subset(clust_df, Dataset == "TCGA")

# format the expr table
res = read_format_MA_expr(tcga_dta, tcga_metadata_table)
tcga_df = res[[1]]

# format table
curr_IDs = tcga_df$ID
gene_ids = colnames(tcga_df)[9:ncol(tcga_df)]
tcga_matr = as.data.frame(t(tcga_df[, 9:ncol(tcga_df)]))
colnames(tcga_matr) = curr_IDs

tcga_matr$hgnc_symbol = gene_ids

# get consensus clusters
tcga_res = get_cluster_labels(gene_map, tcga_matr, tcga_metadata_table, min_consensus_cutoff)

table(tcga_res$verhaak_clust, tcga_res$ClusterK4_kmeans)
```

    ##                
    ##                   1   2   3   4
    ##   IMR_consensus  19   2  97  21
    ##   DIF_consensus   3   9  16 118
    ##   PRO_consensus   1  91   0   8
    ##   MES_consensus 102   0   3   9

``` r
tcga_res$mapped = unlist(lapply(tcga_res$verhaak_clust, mapping_tcga_way))
tcga_consensus_way = tcga_res$ID[which(tcga_res$mapped == tcga_res$ClusterK4_kmeans)]

# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(tcga_res$mapped), 
                                reference = as.factor(tcga_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   1   2   3   4
    ##          1 102   0   3   9
    ##          2   1  91   0   8
    ##          3  19   2  97  21
    ##          4   3   9  16 118
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.8176          
    ##                  95% CI : (0.7809, 0.8506)
    ##     No Information Rate : 0.3126          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.7554          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.005369        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.8160   0.8922   0.8362   0.7564
    ## Specificity            0.9679   0.9773   0.8903   0.9184
    ## Pos Pred Value         0.8947   0.9100   0.6978   0.8082
    ## Neg Pred Value         0.9403   0.9724   0.9472   0.8924
    ## Prevalence             0.2505   0.2044   0.2325   0.3126
    ## Detection Rate         0.2044   0.1824   0.1944   0.2365
    ## Detection Prevalence   0.2285   0.2004   0.2786   0.2926
    ## Balanced Accuracy      0.8920   0.9347   0.8633   0.8374

``` r
tcga_ba = conf_res$byClass[,11]
```

## run mayo

``` r
mayo_data = load(file=file.path(proj_dir, 
                                "/data/mayo/MayoEset.Rda"))

ExpressionData <- get(mayo_data)
mayo_dta <- exprs(ExpressionData)

# read in the metadata
mayo_metadata_table = subset(clust_df, Dataset == "mayo.eset")

# format the expr table
res = read_format_MA_expr(mayo_dta, mayo_metadata_table)
mayo_df = res[[1]]


# format table
curr_IDs = mayo_df$ID
gene_ids = colnames(mayo_df)[9:ncol(mayo_df)]
mayo_matr = as.data.frame(t(mayo_df[, 9:ncol(mayo_df)]))
colnames(mayo_matr) = curr_IDs

mayo_matr$hgnc_symbol = gene_ids


# get consensus clusters
mayo_res = get_cluster_labels(gene_map, mayo_matr, mayo_metadata_table, min_consensus_cutoff)

table(mayo_res$verhaak_clust, mayo_res$ClusterK4_kmeans)
```

    ##                
    ##                  1  2  3  4
    ##   IMR_consensus 11  2 87 17
    ##   DIF_consensus  8 10  1 74
    ##   PRO_consensus  3 65  0  2
    ##   MES_consensus 83  2  5  7

``` r
mayo_res$mapped = unlist(lapply(mayo_res$verhaak_clust, mapping_tcga_way))
mayo_consensus_way = mayo_res$ID[which(mayo_res$mapped == mayo_res$ClusterK4_kmeans)]

# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(mayo_res$mapped), 
                                reference = as.factor(mayo_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 83  2  5  7
    ##          2  3 65  0  2
    ##          3 11  2 87 17
    ##          4  8 10  1 74
    ## 
    ## Overall Statistics
    ##                                          
    ##                Accuracy : 0.8196         
    ##                  95% CI : (0.777, 0.8571)
    ##     No Information Rate : 0.2785         
    ##     P-Value [Acc > NIR] : < 2.2e-16      
    ##                                          
    ##                   Kappa : 0.7587         
    ##                                          
    ##  Mcnemar's Test P-Value : 0.0005065      
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.7905   0.8228   0.9355   0.7400
    ## Specificity            0.9485   0.9832   0.8944   0.9314
    ## Pos Pred Value         0.8557   0.9286   0.7436   0.7957
    ## Neg Pred Value         0.9214   0.9544   0.9769   0.9085
    ## Prevalence             0.2785   0.2095   0.2467   0.2653
    ## Detection Rate         0.2202   0.1724   0.2308   0.1963
    ## Detection Prevalence   0.2573   0.1857   0.3103   0.2467
    ## Balanced Accuracy      0.8695   0.9030   0.9149   0.8357

``` r
mayo_ba = conf_res$byClass[,11]
```

## run tothill

``` r
tothill_data = data(GSE9891_eset)
ExpressionData <- get(tothill_data)
tothill_dta <- exprs(ExpressionData)

tothill_metadata_table = subset(clust_df, Dataset == "Tothill")

# format the expr table
res = read_format_MA_expr(tothill_dta, tothill_metadata_table)
tothill_df = res[[1]]


# format table
curr_IDs = tothill_df$ID
gene_ids = colnames(tothill_df)[9:ncol(tothill_df)]
tothill_matr = as.data.frame(t(tothill_df[, 9:ncol(tothill_df)]))
colnames(tothill_matr) = curr_IDs

tothill_matr$hgnc_symbol = gene_ids

# get consensus clusters
tothill_res = get_cluster_labels(gene_map, tothill_matr, tothill_metadata_table, min_consensus_cutoff)


table(tothill_res$verhaak_clust, tothill_res$ClusterK4_kmeans)
```

    ##                
    ##                  1  2  3  4
    ##   IMR_consensus 10  2 48  5
    ##   DIF_consensus  5  5 11 49
    ##   PRO_consensus  4 38  3  0
    ##   MES_consensus 56  0  5  0

``` r
tothill_res$mapped = unlist(lapply(tothill_res$verhaak_clust, mapping_tcga_way))
tothill_consensus_way = tothill_res$ID[which(tothill_res$mapped == tothill_res$ClusterK4_kmeans)]

# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(tothill_res$mapped), 
                                reference = as.factor(tothill_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 56  0  5  0
    ##          2  4 38  3  0
    ##          3 10  2 48  5
    ##          4  5  5 11 49
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7925          
    ##                  95% CI : (0.7358, 0.8419)
    ##     No Information Rate : 0.3112          
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.722           
    ##                                           
    ##  Mcnemar's Test P-Value : 0.005947        
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.7467   0.8444   0.7164   0.9074
    ## Specificity            0.9699   0.9643   0.9023   0.8877
    ## Pos Pred Value         0.9180   0.8444   0.7385   0.7000
    ## Neg Pred Value         0.8944   0.9643   0.8920   0.9708
    ## Prevalence             0.3112   0.1867   0.2780   0.2241
    ## Detection Rate         0.2324   0.1577   0.1992   0.2033
    ## Detection Prevalence   0.2531   0.1867   0.2697   0.2905
    ## Balanced Accuracy      0.8583   0.9044   0.8094   0.8976

``` r
tothill_ba = conf_res$byClass[,11]
```

## run yoshihara

``` r
yoshihara_data = data(GSE32062.GPL6480_eset)
ExpressionData <- get(yoshihara_data)
yoshihara_dta <- exprs(ExpressionData)

yoshihara_metadata_table = subset(clust_df, Dataset == "Yoshihara")

# format the expr table
res = read_format_MA_expr(yoshihara_dta, yoshihara_metadata_table)
yoshihara_df = res[[1]]


# format table
curr_IDs = yoshihara_df$ID
gene_ids = colnames(yoshihara_df)[9:ncol(yoshihara_df)]
yoshihara_matr = as.data.frame(t(yoshihara_df[, 9:ncol(yoshihara_df)]))
colnames(yoshihara_matr) = curr_IDs

yoshihara_matr$hgnc_symbol = gene_ids


# get consensus clusters
yoshihara_res = get_cluster_labels(gene_map, yoshihara_matr, yoshihara_metadata_table, min_consensus_cutoff)


table(yoshihara_res$verhaak_clust, yoshihara_res$ClusterK4_kmeans)
```

    ##                
    ##                  1  2  3  4
    ##   IMR_consensus 20  0 41  2
    ##   DIF_consensus  3  1 27 57
    ##   PRO_consensus  1 26  3  4
    ##   MES_consensus 65  3  0  2

``` r
yoshihara_res$mapped = unlist(lapply(yoshihara_res$verhaak_clust, mapping_tcga_way))
yoshihara_consensus_way = yoshihara_res$ID[which(yoshihara_res$mapped == yoshihara_res$ClusterK4_kmeans)]

all_consensus = c(schB_consensus_way, schW_consensus_way, tcga_consensus_way,
                  mayo_consensus_way, tothill_consensus_way, yoshihara_consensus_way)

# get the comparison 
conf_res <- caret::confusionMatrix(data=as.factor(yoshihara_res$mapped), 
                                reference = as.factor(yoshihara_res$ClusterK4_kmeans))
conf_res
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction  1  2  3  4
    ##          1 65  3  0  2
    ##          2  1 26  3  4
    ##          3 20  0 41  2
    ##          4  3  1 27 57
    ## 
    ## Overall Statistics
    ##                                           
    ##                Accuracy : 0.7412          
    ##                  95% CI : (0.6828, 0.7938)
    ##     No Information Rate : 0.349           
    ##     P-Value [Acc > NIR] : < 2.2e-16       
    ##                                           
    ##                   Kappa : 0.6463          
    ##                                           
    ##  Mcnemar's Test P-Value : 1.452e-08       
    ## 
    ## Statistics by Class:
    ## 
    ##                      Class: 1 Class: 2 Class: 3 Class: 4
    ## Sensitivity            0.7303   0.8667   0.5775   0.8769
    ## Specificity            0.9699   0.9644   0.8804   0.8368
    ## Pos Pred Value         0.9286   0.7647   0.6508   0.6477
    ## Neg Pred Value         0.8703   0.9819   0.8438   0.9521
    ## Prevalence             0.3490   0.1176   0.2784   0.2549
    ## Detection Rate         0.2549   0.1020   0.1608   0.2235
    ## Detection Prevalence   0.2745   0.1333   0.2471   0.3451
    ## Balanced Accuracy      0.8501   0.9156   0.7289   0.8569

``` r
yoshihara_ba = conf_res$byClass[,11]
```

## format table

``` r
all_dfs = list("SchildkrautB"=schB_ba, "SchildkrautW"=schW_ba, "TCGA"=tcga_ba, "Mayo"=mayo_ba, "Tothill"=tothill_ba, "Yoshihara"=yoshihara_ba)
all_dfs = bind_rows(all_dfs, .id = "Dataset")

colnames(all_dfs) = c("Dataset", "Mesenchymal", "Proliferative", "Immunoreactive", "Differentiated")

library(knitr)
kable(all_dfs, 1, align = "c", digits=3, caption="Balanced Accuracy between consensusOV and our subtype classifications")
```

|   Dataset    | Mesenchymal | Proliferative | Immunoreactive | Differentiated |
|:------------:|:-----------:|:-------------:|:--------------:|:--------------:|
| SchildkrautB |    0.727    |     0.763     |     0.710      |     0.621      |
| SchildkrautW |    0.682    |     0.818     |     0.589      |     0.701      |
|     TCGA     |    0.892    |     0.935     |     0.863      |     0.837      |
|     Mayo     |    0.870    |     0.903     |     0.915      |     0.836      |
|   Tothill    |    0.858    |     0.904     |     0.809      |     0.898      |
|  Yoshihara   |    0.850    |     0.916     |     0.729      |     0.857      |

Balanced Accuracy between consensusOV and our subtype classifications
