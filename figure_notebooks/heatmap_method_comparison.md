recreate_figs_1\_2
================
Natalie Davidson
8/15/2021

# Recreate Figure 1 from Way et.al 2016 –AACES only

Figure 1 from the Way paper shows consistency of each cluster across
datasets. We do the same, but focus on the AACES dataset. This is for
easier visualization for the main figures. We do this same analysis for
K=2, 3, 4.

## First get file paths

``` r
fig1_dir = file.path(proj_dir, 
                    "/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-AcrossCor/")

fig1_files = list.files(fig1_dir,full.names = T )

# only doing this for kmeans, so remove nmf files
fig1_files = grep("nmf", fig1_files, value=T, invert=T)

# lets do the commongenes instead of MAD genes
fig1_files = grep("commongenes", fig1_files, value=T)
```

## Plotting method

``` r
map_samples <- function(in_name){
    
    name_conv = data.frame(in_name = c("aaces.rnaseq","aaces.white.rnaseq", "TCGA", "mayo", "GSE32062.GPL6480", "GSE9891"),
                           out_name = c("Schildkraut B", "Schildkraut W", "TCGA", "Mayo", "Yoshihara", "Tothill"))
    return(name_conv$out_name[which(name_conv$in_name == in_name)]) 
}


plot_single_matr <- function(curr_file, samp_interest){
        
    # the matrix is sample x sample, so len(sample)^2
    list_plots = vector('list', length(samp_interest)^2)
    
    # read in the file
    in_df = data.frame(fread(curr_file))
    row.names(in_df) = in_df$V1
    in_df = in_df[,-1]
    
    in_idx = 1
    for(samp_idx in 2:length(samp_interest)){
        
        samp1 = samp_interest[1]
        samp2 = samp_interest[samp_idx]
    
        
        # get the comparisons of interest
        col_idx = grep(samp1, colnames(in_df))
        row_idx = grep(samp2, rownames(in_df))
        
        curr_df = in_df[row_idx, col_idx]
        colnames(curr_df) = 1:ncol(curr_df)
        curr_df$compare_samp = row.names(curr_df) 

        
        # plot the correlation matrix as a heatmap
        gg_out = ggplot(data = melt(curr_df)) + geom_tile(aes(y=variable,x=compare_samp, fill = value)) +
                        geom_text(aes(x=compare_samp,y=variable, label=round(value, 2))) + 
                        scale_fill_gradient2(low="blue", mid="white", high="red", 
                                            midpoint=0,    
                                            breaks=seq(-1,1,0.1), 
                                            limits=c(-1, 1)) + 
                        theme_bw() + labs(y = map_samples(samp1), x=map_samples(samp2)) +
                        theme(legend.position = "none", 
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.text.y=element_blank(),
                                axis.ticks.y=element_blank())
        # since we make several plot stacked along a horizontal line, 
        # only keep the y-axis title of the first one
        if(in_idx != 1){
            gg_out = gg_out + theme(axis.title.y=element_blank())
        }
        

        list_plots[[in_idx]] = gg_out
        in_idx = in_idx +1 

    }

    flat_plot = do.call(ggarrange, 
                            c(list_plots, 
                                ncol=length(samp_interest)-1,
                                nrow=1))

    return(flat_plot)
}
```

## Now for each file make the corr. plot

``` r
# we are going to sort the filenames
# this way we can be sure the results
# are in the following order K=2, 3, 4
fig1_files = sort(fig1_files)

# we will only plot the sample comparisons of interest
samp_interest = c("aaces.rnaseq","aaces.white.rnaseq", "TCGA", "mayo", "GSE32062.GPL6480", "GSE9891")
k3_matr = plot_single_matr(fig1_files[2], samp_interest)


outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K3_method_comparison_AA.pdf")
ggsave(outfile,
       k3_matr[[1]], width = 15, height = 3, units = "in", device = "pdf")


# now re run again for K=4
k4_matr = plot_single_matr(fig1_files[3], samp_interest)


outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K4_method_comparison_AA.pdf")
ggsave(outfile,
       k4_matr[[1]], width = 15, height = 3, units = "in", device = "pdf")



# we will only plot the sample comparisons of interest
samp_interest = c("aaces.white.rnaseq","aaces.rnaseq", "TCGA", "mayo", "GSE32062.GPL6480", "GSE9891")
k3_matr = plot_single_matr(fig1_files[2], samp_interest)


outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K3_method_comparison_W.pdf")
ggsave(outfile,
       k3_matr[[1]], width = 15, height = 3, units = "in", device = "pdf")


# now re run again for K=4
k4_matr = plot_single_matr(fig1_files[3], samp_interest)


outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K4_method_comparison_W.pdf")
ggsave(outfile,
       k4_matr[[1]], width = 15, height = 3, units = "in", device = "pdf")
```
