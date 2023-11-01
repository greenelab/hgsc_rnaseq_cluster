recreate_figs_1\_2
================
Natalie Davidson
8/15/2021

# Recreate Figure 1 from Way et.al 2016

Figure 1 from the Way paper shows consistency of each cluster across
datasets. We do this same analysis for K=2, 3, 4.

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
    
    name_conv = data.frame(in_name = c("aaces.rnaseq", "aaces.white.rnaseq", "TCGA", "mayo", "GSE32062.GPL6480", "GSE9891"),
                           out_name = c("Schildtkraut B","Schildtkraut W",  "TCGA", "Mayo", "Yoshihara", "Tothill"))
    return(name_conv$out_name[which(name_conv$in_name == in_name)]) 
}



plot_triangle_matr <- function(curr_file, samp_interest){
        
    # the matrix is sample x sample, so len(sample)^2
    list_plots = vector('list', length(samp_interest)^2)
    
    # read in the file
    in_df = data.frame(fread(curr_file))
    row.names(in_df) = in_df$V1
    in_df = in_df[,-1]
    
    in_idx = 1
    for(samp1_idx in 1:length(samp_interest)){
        for(samp2_idx in 1:length(samp_interest)){
            
            # this check enforces lower triangular
            if(samp1_idx >= samp2_idx){
                in_idx = in_idx +1 
                next
            }else{
                samp1 = samp_interest[samp1_idx]
                samp2 = samp_interest[samp2_idx]
            }
            
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
                            theme_bw() + labs(x = map_samples(samp2), y=map_samples(samp1)) +
                            theme(legend.position = "none", 
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.ticks.y=element_blank())
            if(samp1_idx != 1){
                gg_out = gg_out + theme(axis.title.x=element_blank())
            }
            if(samp2_idx != 6){
                gg_out = gg_out + theme(axis.title.y=element_blank())
            }
            list_plots[[in_idx]] = gg_out
            in_idx = in_idx +1 
    
        }
    }
    triangle_plot = do.call(ggarrange, 
                            c(rev(list_plots), 
                                nrow=length(samp_interest),
                                ncol=length(samp_interest)))

    return(triangle_plot)
}
```

## Now for each file make the corr. plot

``` r
# we are going to sort the filenames
# this way we can be sure the results
# are in the following order K=2, 3, 4
fig1_files = sort(fig1_files)

# we will only plot the sample comparisons of interest
samp_interest = rev(c("aaces.rnaseq", "aaces.white.rnaseq", "TCGA", "mayo", "GSE32062.GPL6480", "GSE9891"))

k2_matr = plot_triangle_matr(fig1_files[1], samp_interest)
k3_matr = plot_triangle_matr(fig1_files[2], samp_interest)
k4_matr = plot_triangle_matr(fig1_files[3], samp_interest)


outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K2_full_method_comparison.pdf")
ggsave(outfile,
       k2_matr, width = 15, height = 15, units = "in", device = "pdf")
outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K3_full_method_comparison.pdf")
ggsave(outfile,
       k3_matr, width = 15, height = 15, units = "in", device = "pdf")
outfile = paste0(proj_dir, "/figure_notebooks/manuscript_figs/K4_full_method_comparison.pdf")
ggsave(outfile,
       k4_matr, width = 15, height = 15, units = "in", device = "pdf")
```
