Tech\_comparison
================
Natalie Davidson
3/19/2021

## Overview

This document will compare the expression profiles across 3 technologies: RNA-Seq, Nanostring, and HTA To do this we will need the expression tables and sample ID mapping tables. The data is located in 2 repos: [hgsc\_characterization](https://github.com/greenelab/hgsc_characterization/) and [aaces\_ovarian](https://github.com/greenelab/aaces_ovarian/)

To re-run this analysis, it is assumed that you are in the hgsc\_characterization Rproj and both repos are in a shared folder.

**Expression Table file locations:**

-   Nanostring expression:
    -   repo: [hgsc\_characterization](https://github.com/greenelab/hgsc_characterization/)
    -   location: `/hgsc_characterization/data/nanostring/NCONanostingClassifierData_forAriel.csv`
    -   notes: not all gene expression is reported. There is also very few samples.
-   RNA-Seq expression:
    -   repo: [hgsc\_characterization](https://github.com/greenelab/hgsc_characterization/)
    -   location: `/hgsc_characterization/data/rna_seq/salmon_quant_processed/salmon_gene_quant.RDS`
-   HTA expression:
    -   repo: [aaces\_ovarian](https://github.com/greenelab/aaces_ovarian/)
    -   location: `/aaces_ovarian/output/AACES_expression_scanned.txt`

**Sample ID mapping**. The sample ID used to index all samples in this analysis is the **suid**.

-   Nanostring expression:
    -   no action needed
-   HTA sample IDs
    -   HTA sample IDs look like: 5716\_aaaaa\_aaaaaa\_HTA.CEL where a is 0-9.
    -   To translate HTA sample ids to **suid**: `/aaces_ovarian/data/exclusion_data.tsv`
-   RNA-Seq expression:
    -   RNA-Seq sample IDs look like 18341Xaa where a is 0-9.
    -   To translate RNA-Seq sample ids to **suid**: `/hgsc_characterization/reference_data/rna_sample_metadata.txt`

**Sample Type Annotation**.

-   `/hgsc_characterization/data/nanostring_hta_overlap_samples.csv`
    -   The sample subtype used to find sample patterns
    -   **HGSCsubtypeFinal** is the column of interest.

**Helper Scripts**. helper scripts + methods are listed here:

-   `/hgsc_characterization/comparison/utils/file_processing_utils.R`
    -   format\_nanostring\_data
    -   format\_hta\_data
    -   format\_rnaseq\_data
    -   get\_gene\_id\_map
-   `/hgsc_characterization/analyze_rnaseq/plot_utils.R`
    -   display\_venn
    -   plot\_pca

## Read in files

Read each data source

``` r
# load expression data
nano_expr_file = file.path(proj_dir, 
                    "/data/nanostring/NCONanostingClassifierData_forAriel.csv")
nano_dt = format_nanostring_data(nano_expr_file)

hta_expr_file = file.path(proj_dir, 
                    "../aaces_ovarian/output/AACES_expression_scanned.txt")
hta_trans_file = file.path(proj_dir,
                           "../aaces_ovarian/data/exclusion_data.tsv")
hta_dt = format_hta_data(hta_expr_file, hta_trans_file)


rnaseq_expr_file = file.path(proj_dir, 
                    "/data/rna_seq/salmon_quant_processed/salmon_gene_quant.RDS")

rnaseq_trans_file = file.path(proj_dir, 
                    "/reference_data/rna_sample_metadata.txt")
rnaseq_dt = format_rnaseq_data(rnaseq_expr_file, rnaseq_trans_file)

# sample subtype annotation
subtype_file = file.path(proj_dir, 
                    "/data/nanostring_hta_overlap_samples.csv")
subtype_dt = fread(subtype_file)
subtype_dt = unique(subtype_dt[,c("suid", "HGSCsubtypeFinal")])
```

## Basic QC

Get the intersection of gene IDs

``` r
nano_hgnc = unique(na.omit(nano_dt$hgnc_symbol))
hta_hgnc = unique(na.omit(hta_dt$hgnc_symbol))
rnaseq_hgnc = unique(na.omit(rnaseq_dt$hgnc_symbol))

gene_per_tech = list(nano = nano_hgnc,
                     hta = hta_hgnc,
                     rnaseq = rnaseq_hgnc)
display_venn(gene_per_tech,
            category.names = c("Nanostring" , "HTA", "RNA-Seq"),
            fill = c("#1B9E77", "#D95F02", "#7570B3"),
            main = "Gene Overlap across 3 Technologies")
```

<img src="tech_comparison_files/figure-markdown_github/get-gene_overlap-1.png" style="display: block; margin: auto;" />

``` r
# do they all have the UBC genes?
ubc_file = file.path(proj_dir, "/reference_data/UBC_genes.txt")
ubc_df <- fread(ubc_file)
colnames(ubc_df)[1] = "hgnc_symbol"
gene_ids_df = get_gene_id_map(unique(ubc_df$hgnc_symbol), 
                                filter_type="hgnc_symbol",
                                attributes= c("hgnc_symbol", "ensembl_gene_id"))
ubc_df = merge(as.data.table(gene_ids_df), ubc_df, all.y=T, by="hgnc_symbol")


gene_per_tech = list(nano = nano_hgnc,
                     hta = hta_hgnc,
                     rnaseq = rnaseq_hgnc,
                     ubc=unique(ubc_df$hgnc_symbol))
display_venn(gene_per_tech,
            category.names = c("Nanostring" , "HTA", "RNA-Seq", "UBC"),
            fill = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
            main = "Gene Overlap")
```

<img src="tech_comparison_files/figure-markdown_github/get-gene_overlap-2.png" style="display: block; margin: auto;" />

``` r
# write out gene names is missing
setdiff(ubc_df$hgnc_symbol, nano_hgnc)
```

    ## character(0)

``` r
setdiff(ubc_df$hgnc_symbol, hta_hgnc)
```

    ## character(0)

``` r
setdiff(ubc_df$hgnc_symbol, rnaseq_hgnc)
```

    ## character(0)

Get the intersection of sample IDs

``` r
nano_suid = unique(na.omit(nano_dt$suid))
hta_suid = unique(na.omit(hta_dt$suid))
rnaseq_suid = unique(na.omit(rnaseq_dt$suid))

samp_per_tech = list(nano = nano_suid,
                     hta = hta_suid,
                     rnaseq = rnaseq_suid)
display_venn(samp_per_tech,
            category.names = c("Nanostring" , "HTA", "RNA-Seq"),
            fill = c("#1B9E77", "#D95F02", "#7570B3"),
            main = "Sample Overlap across 3 Technologies")
```

<img src="tech_comparison_files/figure-markdown_github/get_samp_overlap-1.png" style="display: block; margin: auto;" />

## Correlations across Technologies

Now lets take the samples and genes that are in all 3 technologies and compare them to one another.

``` r
genes_in_all = Reduce(intersect, gene_per_tech[1:3])
samps_in_all = Reduce(intersect, samp_per_tech)

nano_dt$tech = "nanostring"
hta_dt$tech = "hta"
rnaseq_dt$tech = "rnaseq"

total_dt = rbind(nano_dt, hta_dt, rnaseq_dt)
total_dt = merge(subtype_dt, total_dt, by="suid")

intersect_dt = subset(total_dt, 
                      suid %in% samps_in_all &
                      hgnc_symbol %in% genes_in_all)

# to make technologies comparable scale within each sample
intersect_dt_scaled = intersect_dt %>%
                        group_by(suid, tech) %>%
                        mutate(scaled_expr = scale(expr))
                        

# now dcast so we can analyze them
intersect_dt_cast = dcast(intersect_dt_scaled, suid+tech+HGSCsubtypeFinal ~ hgnc_symbol, value.var="scaled_expr")
row.names(intersect_dt_cast) = paste(intersect_dt_cast$suid,
                                     intersect_dt_cast$tech, sep="_")

# now plot heatmap across all genes
expr_df = subset(intersect_dt_cast, select=-c(suid, tech, HGSCsubtypeFinal))

annotation_row_df <- data.frame(subset(intersect_dt_cast, select=c(tech, HGSCsubtypeFinal)))
annotation_row_df$HGSCsubtypeFinal = as.factor(annotation_row_df$HGSCsubtypeFinal)
pheatmap(expr_df, annotation_row = annotation_row_df, 
         annotation_col = annotation_row_df, 
         scale = "none", show_rownames = T,
         show_colnames = F, cluster_cols = T, cluster_rows = T)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps-1.png" style="display: block; margin: auto;" />

``` r
# now plot heatmap of correlation matrix
cor_df <- cor(t(expr_df), method="spearman")
pheatmap(cor_df, annotation_row = annotation_row_df, 
         annotation_col = annotation_row_df, 
         scale = "none", show_rownames = T,
         show_colnames = F, cluster_cols = T, cluster_rows = T)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps-2.png" style="display: block; margin: auto;" />

``` r
# make a line plot to do pairwise comparisons
cor_line_df = cbind(subset(intersect_dt_cast, select=c(suid, tech)),
                    cor_df)
cor_line_df = melt(cor_line_df, id.vars = c("suid", "tech"))
cor_line_df = cor_line_df[ grep("_rnaseq", cor_line_df$variable ), ]
cor_line_df = subset(cor_line_df, tech != "rnaseq")
colnames(cor_line_df)[3:4] = c("cor_id", "spearman_cor")

ggplot(cor_line_df, aes(x=as.factor(suid), y=spearman_cor, fill=tech)) +
    geom_bar(stat = "identity", position="dodge") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
    xlab("Sample IDs") + ylab("Spearman Correlation with RNA-Seq") +
    ggtitle("Spearman correlation with RNA-Seq, compared to HTA and Nanostring (Higher is better)") +
    scale_fill_brewer(palette="Set2")
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps-3.png" style="display: block; margin: auto;" />

``` r
# now a PCA to see if everything is far off
plot_pca(t(expr_df), title = "PCA across technologies", scale=F)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps-4.png" style="display: block; margin: auto;" />

## Correlations across Technologies Restricted to UBC genes

Now lets restrict to UBC genes to see if it is any better.

``` r
# this now includes an intersection with UBC genes
genes_in_all = Reduce(intersect, gene_per_tech)
intersect_dt = subset(total_dt, 
                        suid %in% samps_in_all &
                        hgnc_symbol %in% genes_in_all)

# to make technologies comparable scale within each sample
intersect_dt_scaled = intersect_dt %>%
                        group_by(suid, tech) %>%
                        mutate(scaled_expr = scale(expr))
                        

# now dcast so we can analyze them
intersect_dt_cast = dcast(intersect_dt_scaled, suid+tech+HGSCsubtypeFinal ~ hgnc_symbol, value.var="scaled_expr")
row.names(intersect_dt_cast) = paste(intersect_dt_cast$suid,
                                     intersect_dt_cast$tech, sep="_")

# now plot heatmap across all genes
expr_df = subset(intersect_dt_cast, select=-c(suid, tech, HGSCsubtypeFinal))

annotation_row_df <- subset(intersect_dt_cast, select=c(tech, HGSCsubtypeFinal))
annotation_row_df$HGSCsubtypeFinal <- as.factor(annotation_row_df$HGSCsubtypeFinal)
pheatmap(expr_df, annotation_row = annotation_row_df, 
         scale = "none", show_rownames = T,
         show_colnames = F, cluster_cols = T, cluster_rows = T)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps_UBC-1.png" style="display: block; margin: auto;" />

``` r
# now plot heatmap of correlation matrix
cor_df <- cor(t(expr_df), method="spearman")
pheatmap(cor_df, annotation_row = annotation_row_df, 
          annotation_col = annotation_row_df, 
         scale = "none", show_rownames = T,
         show_colnames = F, cluster_cols = T, cluster_rows = T)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps_UBC-2.png" style="display: block; margin: auto;" />

``` r
# make a line plot to do pairwise comparisons
cor_line_df = cbind(subset(intersect_dt_cast, select=c(suid, tech)),
                    cor_df)
cor_line_df = melt(cor_line_df, id.vars = c("suid", "tech"))
cor_line_df = cor_line_df[ grep("_rnaseq", cor_line_df$variable ), ]
cor_line_df = subset(cor_line_df, tech != "rnaseq")
colnames(cor_line_df)[3:4] = c("cor_id", "spearman_cor")

ggplot(cor_line_df, aes(x=as.factor(suid), y=spearman_cor, fill=tech)) +
    geom_bar(stat = "identity", position="dodge") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
    xlab("Sample IDs") + ylab("Spearman Correlation with RNA-Seq") +
    ggtitle("Spearman correlation with RNA-Seq, compared to HTA and Nanostring (Higher is better)") +
    scale_fill_brewer(palette="Set2")
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps_UBC-3.png" style="display: block; margin: auto;" />

``` r
# now a PCA to see if everything is far off
plot_pca(t(expr_df), title = "PCA across technologies", scale=F)
```

<img src="tech_comparison_files/figure-markdown_github/get_orrelation_heatmaps_UBC-4.png" style="display: block; margin: auto;" />
