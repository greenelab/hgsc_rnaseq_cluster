library("ggplot2")

cluster_purity <- function(clusters, classes) {
  # taken from https://stackoverflow.com/questions/9253843/r-clustering-purity-metric
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}


plot_pca <- function(in_df, title, color_df=NA, color_id="Subtype", scale=TRUE){

    pca_df = t(data.frame(in_df))
    keep_idx = apply(pca_df, 2, function(x) var(x) != 0)
    pca_df = pca_df[,keep_idx]

    pca_res = prcomp(pca_df, scale. = scale)
    pca_res_df = data.frame(pca_res$x)
    pca_res_df$Sample = row.names(pca_res_df)

    percentage <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 2)
    percentage <- paste( colnames(pca_res_df), "(", paste( as.character(percentage), "%", ")", sep="") )

    # color by sample if asked for
    if(!is.na(color_df)){

        if(length(grep("Sample", colnames(color_df))) == 0){
            stop("to use the color ability, you need a column id called Sample")
        }

        pca_res_df = merge(color_df, pca_res_df, by="Sample")
        gg_pca = ggplot(pca_res_df, aes_string(x="PC1",y="PC2", color=color_id)) +
            geom_point()  +
            theme_bw() + xlab(percentage[1]) + ylab(percentage[2]) +
            ggtitle(title)
    }else{

        gg_pca = ggplot(pca_res_df,
                    aes_string(x="PC1",y="PC2",label="Sample")) +
            geom_point()  +
            theme_bw() + xlab(percentage[1]) + ylab(percentage[2]) +
            ggtitle(title) + geom_text_repel()
    }




    return(gg_pca)
}

plot_umap <- function(in_df, title, color_df=NA, color_id="Subtype", scale=TRUE){

    umap_df = t(data.frame(in_df))
    keep_idx = apply(umap_df, 2, function(x) var(x) != 0)
    umap_df = umap_df[,keep_idx]

    umap_res = umap(umap_df)
    umap_res_df = data.frame(umap_res$layout)
    colnames(umap_res_df) = c("umap1", "umap2")
    umap_res_df$Sample = row.names(umap_res_df)

    # color by sample if asked for
    if(!is.na(color_df)){

        if(length(grep("Sample", colnames(color_df))) == 0){
            stop("to use the color ability, you need a column id called Sample")
        }

        umap_res_df = merge(color_df, umap_res_df, by="Sample")
        gg_umap = ggplot(umap_res_df,
                    aes_string(x="umap1",y="umap2",
                               label="Sample", color=color_id)) +
            theme_bw() + xlab("UMAP 1") + ylab("UMAP 2") +
            ggtitle(title) + geom_text_repel()
    }else{

        gg_umap = gg_umap+geom_point()  +
            theme_bw() + xlab("UMAP 1") + ylab("UMAP 2") +
            ggtitle(title) + geom_text_repel()
    }




    return(gg_umap)
}


display_venn <- function(x, ...){
  library("VennDiagram")
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#' Method to plot a confusion matrix
#' for a multi-class prediction task
#' Actual and Predict must be in the same order
#' and they must have the same number of levels
#'
#' @param Actual, true labels
#' @param Predict, predicted labels
#' @param colors, color scale palette (mid, low, high)
#' @param text.scl, scaled size of text in the confusion matrix
#' @return ggplot object
prettyConfused <- function(Actual, Predict, colors=c("white","red4","dodgerblue3"), text.scl=5){

    #### taken from https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot
    actual = as.data.frame(table(Actual))
    names(actual) = c("Actual","ActualFreq")

    #build confusion matrix
    confusion = as.data.frame(table(Actual, Predict))
    names(confusion) = c("Actual","Predicted","Freq")

    #calculate percentage of test cases based on actual frequency

    confusion = merge(confusion, actual, by=c('Actual','Actual'))
    confusion$Percent = confusion$Freq/(confusion$ActualFreq+0.000001)*100
    confusion$ColorScale<-confusion$Percent*-1
    confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
    confusion$Label<-paste(round(confusion$Percent,0),"%, \nn=",confusion$Freq,sep="")
    tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")

    tile = tile +
        geom_text(aes(x=Actual,y=Predicted, label=Label),data=confusion, size=text.scl, colour="black") +
        scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0,guide='none')

    return(tile)
}
