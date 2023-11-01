
######################################
### In this script we calculate
### the RNA degradation of Agilent microarray samples
### this is done using the arrayQuality
### package, provides basic QC measures
######################################

detach("package:ggpubr", unload = TRUE)
detach("package:rstatix", unload = TRUE)
detach("package:broom", unload = TRUE)
detach("package:tidyr", unload = TRUE)
detach("package:dplyr", unload = TRUE)
detach("package:cowplot", unload = TRUE)
detach("package:ggplot2", unload = TRUE)

suppressPackageStartupMessages({
    library(marray)
    library(arrayQuality)
})


yoshihara_files = "/Users/davidnat/Downloads/GSE32062_RAW/"
yoshihara_files = list.files(yoshihara_files, full.names = F)
yoshihara_files = grep(".txt.gz", yoshihara_files, value=T, invert=F)
setwd("/Users/davidnat/Downloads/GSE32062_RAW/")
# I think this is single channel, so just duplicate for green and red so they are the same
# idea taken from https://stat.ethz.ch/pipermail/bioconductor/2008-April/022052.html
agil_data = read.Agilent(fnames = yoshihara_files[1],
                         path=NULL,
                         name.Gf = "gMedianSignal", name.Gb = "gBGMedianSignal",
                         name.Rf = "gMedianSignal", name.Rb = "gBGMedianSignal",
                         name.W= NULL, layout = NULL, gnames = NULL, targets = NULL,
                         notes=NULL, skip=NULL, sep="\t", quote="\"", DEBUG=FALSE,
                         info.id=NULL)
maQualityPlots(agil_data)




mayo_files = "/Users/davidnat/Downloads/GSE74357_RAW/"
mayo_files = list.files(mayo_files, full.names = F)
mayo_files = grep(".txt.gz", mayo_files, value=T, invert=F)
setwd("/Users/davidnat/Downloads/GSE74357_RAW/")
# I think this is single channel, so just duplicate for green and red so they are the same
# idea taken from https://stat.ethz.ch/pipermail/bioconductor/2008-April/022052.html
agil_data = read.Agilent(fnames = mayo_files[1],
                         path=NULL,
                         name.Gf = "gMedianSignal", name.Gb = "gBGMedianSignal",
                         name.Rf = "rMedianSignal", name.Rb = "rBGMedianSignal",
                         name.W= NULL, layout = NULL, gnames = NULL, targets = NULL,
                         notes=NULL, skip=NULL, sep="\t", quote="\"", DEBUG=FALSE,
                         info.id=NULL)
maQualityPlots(agil_data)
qcScore(agil_data)
