######################################
### In this script we calculate
### the RNA degradation of Affy microarray samples
### this is done using the AffyRNADegradation
### package, which calculated the 5'/3' bias
######################################


suppressPackageStartupMessages({
    library("ggpubr")
    library("rstatix")
    library("broom")
    library("tidyr")
    library("dplyr")
    library("cowplot")
    library("ggplot2")
})



detach("package:ggpubr", unload = TRUE)
detach("package:rstatix", unload = TRUE)
detach("package:broom", unload = TRUE)
detach("package:tidyr", unload = TRUE)
detach("package:dplyr", unload = TRUE)
detach("package:cowplot", unload = TRUE)
detach("package:ggplot2", unload = TRUE)

suppressPackageStartupMessages({
    library("AffyRNADegradation")
    library("affy")
})


affy_tcga = affy::ReadAffy(celfile.path="/Users/davidnat/Downloads/GSE26712_RAW/", compress=TRUE)

tongs <- AffyRNADegradation::GetTongs(affy_tcga, chip.idx = 1)
AffyRNADegradation::PlotTongs(tongs)

all_d = c()
all_slope = c()
all_3 = c()
all_5 = c()
all_sampname = c()
for(curr_idx in 1:length(affy_tcga$sample)){
    print(curr_idx)
    startTime <- Sys.time()


    rna.deg <- affy::AffyRNAdeg(affy_tcga[,curr_idx])

    all_slope = c(all_slope, rna.deg$slope)
    all_3 = c(all_3, rna.deg$means.by.number[1])
    all_5 = c(all_5, rna.deg$means.by.number[11])

    curr_name = rna.deg$sample.names
    curr_name = strsplit(curr_name, "_")[[1]][1]
    all_sampname = c(all_sampname, curr_name)
    gc()

    endTime <- Sys.time()
    print(endTime - startTime)
}

slope_df = data.frame(ID=all_sampname,
                      slope=all_slope,
                      threeP=all_3,
                      fiveP=all_5)

bonome_file = file.path("~/Documents/projects/greenelab/hgsc_characterization/",
                         "/reference_data/bonome_slope.tsv")
write.table(slope_df, bonome_file, quote=F, row.names = F, sep="\t")
