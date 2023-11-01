#This script brings in the RNAseq and epidemiological data for 
#all six datasets and harmonizes variables across datasets. 
#Output from this script includes a dataset to be used for 
#survival analyses. 

install.packages("here")
install.packages("survival")
install.packages("rlang")
install.packages("dplyr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedOvarianData")

library("here")
library("survival")
library("rlang")
library("dplyr")
library("curatedOvarianData")

proj_dir = here()

#readin cluster assignments from way pipeline output
clusterinfo = file.path(proj_dir, "data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv")
keepall <-read.csv(clusterinfo, header=TRUE, sep="," )

#Add AACES 
  epidata = file.path(proj_dir, "reference_data/tissuegrant_epidata_07212023.csv")
  aacesrawepi <- read.csv(epidata, header=TRUE, sep=",")
  
  #format epi variables 
  aacesrawepi$age <- as.numeric(aacesrawepi$refage)
  aacesrawepi$cd3 <- as.numeric(aacesrawepi$meancd3)
  aacesrawepi$cd3cd8 <- as.numeric(aacesrawepi$meancd3cd8)
  aacesrawepi$stage <- as.factor(aacesrawepi$figo_stage)
  aacesrawepi$T1stage <- "Missing"
  aacesrawepi$T1stage[aacesrawepi$stage=="I"] <- "I"
  aacesrawepi$T1stage[aacesrawepi$stage=="II"] <- "II"
  aacesrawepi$T1stage[aacesrawepi$stage=="III"] <- "III"
  aacesrawepi$T1stage[aacesrawepi$stage=="IV"] <- "IV"
  
  aacesrawepi$debulking <- as.factor(aacesrawepi$dblkstat)
  levels(aacesrawepi$debulking)[levels(aacesrawepi$debulking)==""] <- "Missing"
  levels(aacesrawepi$debulking)[levels(aacesrawepi$debulking)=="no surgery"] <- "Suboptimal"
  levels(aacesrawepi$debulking)[levels(aacesrawepi$debulking)=="optimal"] <- "Optimal"
  levels(aacesrawepi$debulking)[levels(aacesrawepi$debulking)=="suboptimal"] <- "Suboptimal"
  aacesrawepi$dataset <- "AA"
  aacesrawepi$dataset <- as.factor(aacesrawepi$dataset)
  aacesrawepi$diagyrcat <- "Missing"
  aacesrawepi$diagyrcat[aacesrawepi$diagyear > 1998 & aacesrawepi$diagyear < 2005] <- "1999-2004"
  aacesrawepi$diagyrcat[aacesrawepi$diagyear > 2004 & aacesrawepi$diagyear < 2011] <- "2005-2010"
  aacesrawepi$diagyrcat[aacesrawepi$diagyear > 2010 & aacesrawepi$diagyear < 2016] <- "2011-2015"
  aacesrawepi$study <- as.factor(aacesrawepi$study)

  levels(aacesrawepi$neoadj)[levels(aacesrawepi$neoadj)==""] <- "Missing"
  levels(aacesrawepi$neoadj)[levels(aacesrawepi$neoadj)=="no"] <- "No"
  levels(aacesrawepi$neoadj)[levels(aacesrawepi$neoadj)=="yes"] <- "Yes"
  
  #rename and recode survival variables to match the other datasets
  aacesrawepi$age_at_initial_pathologic_diagnosis <- aacesrawepi$age
  aacesrawepi$days_to_death <- aacesrawepi$survival_days
  aacesrawepi$vital_status[aacesrawepi$vitalstatus=="Dead"] <- "deceased"
  aacesrawepi$vital_status[aacesrawepi$vitalstatus=="Alive"] <- "living"
  
  #link epi data and cluster info
  aametatable = file.path(proj_dir, "reference_data/main_AA_metadata_table.tsv")
  whitemetatable = file.path(proj_dir, "reference_data/main_white_metadata_table.tsv")
  aametaraw <-read.table(aametatable, header=TRUE, sep="\t")
  aametaraw$suid <- substr(aametaraw$suid,1,6)
  whitemetaraw <-read.table(whitemetatable, header=TRUE, sep="\t")
  whitemetaraw$suid <- substr(whitemetaraw$suid,1,5)
  aameta<-aametaraw[,c("suid", "ran_in_way_pipeline", "ClusterK2_kmeans", "ClusterK3_kmeans","ClusterK4_kmeans", "ClusterK4_kmeans_TCGA_names")]
  whitemeta <-whitemetaraw[,c("suid", "ran_in_way_pipeline","ClusterK2_kmeans", "ClusterK3_kmeans","ClusterK4_kmeans", "ClusterK4_kmeans_TCGA_names")]
  aaepi <- merge(aacesrawepi,aameta,by.x="suid", by.y="suid", all.y=TRUE)
  whiteepi <- merge(aacesrawepi,whitemeta,by.x="suid", by.y="suid", all.y=TRUE)

     #remove cases with low expression or otherwise not run in way pipeline
     SchildkrautBlack <- filter(aaepi, aaepi$ran_in_way_pipeline==TRUE)
     SchildkrautWhite <- filter(whiteepi, whiteepi$ran_in_way_pipeline==TRUE)
     #remove cases with neoadjuvant chemo 
     SchildkrautBlack <- filter(SchildkrautBlack, !SchildkrautBlack$neoadj=="Yes")
     SchildkrautWhite <- filter(SchildkrautWhite, !SchildkrautWhite$neoadj=="Yes")
     #keep only black women in SchildkrautBlack and only white women in SchildkrautWhite
     SchildkrautBlack <- filter(SchildkrautBlack, SchildkrautBlack$race=="black")
     SchildkrautWhite <- filter(SchildkrautWhite, SchildkrautWhite$race=="white")
     

  SchildkrautBlack$sampleid <- SchildkrautBlack$suid
  SchildkrautBlack$Dataset <- "SchildkrautB"
  SurvAA <- select(SchildkrautBlack, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)

    TableOneBlack <- select(SchildkrautBlack, sampleid, ClusterK3_kmeans, ClusterK4_kmeans, age, T1stage, debulking, neoadj, cd3, cd3cd8)
    
  SchildkrautWhite$sampleid <- SchildkrautWhite$suid
  SchildkrautWhite$Dataset <- "SchildkrautW"
  SurvWhite <- select(SchildkrautWhite, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)
    
    TableOneWhite <- select(SchildkrautWhite, sampleid, ClusterK3_kmeans, ClusterK4_kmeans, age, T1stage, debulking, neoadj, cd3, cd3cd8)
    
    
#Add Mayo
  epimayo = file.path(proj_dir, "data/mayo/Mayo_Pheno_Data.csv")
  mayorawepi <- read.csv(epimayo, header=TRUE, sep=",")
  #rename and recode survival variables to match the other datasets
    mayorawepi$stage <- as.factor(mayorawepi$tumorstage)
    mayorawepi$debulking <- as.factor(mayorawepi$debulking)
      levels(mayorawepi$debulking)[levels(mayorawepi$debulking)=="optimal"] <- "Optimal"
      levels(mayorawepi$debulking)[levels(mayorawepi$debulking)=="suboptimal"] <- "Suboptimal"
    mayorawepi$sampleid <- mayorawepi$unique_patient_ID
    Survmayo <- merge(mayorawepi, keepall, by.x = "sampleid", by.y="X", all.x=FALSE, all.y=FALSE)
    Survmayo <- select(Survmayo, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)
    Survmayo$Dataset <- "Mayo"


#Get curateOvarianData from bioconductor (TCGA, Toth, Yoshihara)

    args <- c("TCGA_eset", "GSE32062.GPL6480_eset", "GSE9891_eset")
    
  #get data in curatedOvarianData
    detailedData <- data(package="curatedOvarianData")[3]
    detailedData.names <- detailedData$results[,3]
    inCuratedOvarianData <- args %in% detailedData.names
    
    argsCurated <- args[inCuratedOvarianData]
    
    #load the phenotype data (dta)
    phenoData <- list()
    for (phenoset in 1:3) {
      dta <- c()
      if (args[phenoset] %in% detailedData.names) {
      data(list = args[phenoset], package = "curatedOvarianData")
      tmp <- get(args[phenoset])
      phenoData[[phenoset]] <- pData(tmp)
      names(phenoData)[phenoset] <- args[phenoset]
    }}
    
#Add TCGA
  phenoTCGA <- as.data.frame(phenoData["TCGA_eset"])
  phenoTCGA <- tibble::rownames_to_column(phenoTCGA, "sampleid")
  SurvTCGA <- merge(phenoTCGA, keepall, by.x = "sampleid", by.y = "X", all.x=FALSE, all.y=FALSE)

  SurvTCGA$stage <- SurvTCGA$TCGA_eset.tumorstage
  SurvTCGA$debulking <- SurvTCGA$TCGA_eset.debulking
    levels(SurvTCGA$debulking)[levels(SurvTCGA$debulking)=="optimal"] <- "Optimal"
    levels(SurvTCGA$debulking)[levels(SurvTCGA$debulking)=="suboptimal"] <- "Suboptimal"
  SurvTCGA$age_at_initial_pathologic_diagnosis <- SurvTCGA$TCGA_eset.age_at_initial_pathologic_diagnosis
  SurvTCGA$vital_status <- SurvTCGA$TCGA_eset.vital_status
  SurvTCGA$days_to_death <- SurvTCGA$TCGA_eset.days_to_death
  SurvTCGA <- select(SurvTCGA, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)

#Add Tothill
  phenoToth <- as.data.frame(phenoData["GSE9891_eset"])
  phenoToth <- tibble::rownames_to_column(phenoToth, "sampleid")
  SurvToth <- merge(phenoToth, keepall, by.x = "sampleid", by.y = "X", all.x=FALSE, all.y=FALSE)
  
  SurvToth$stage <- SurvToth$GSE9891_eset.tumorstage
  SurvToth$debulking <- SurvToth$GSE9891_eset.debulking
  levels(SurvToth$debulking)[levels(SurvToth$debulking)=="optimal"] <- "Optimal"
  levels(SurvToth$debulking)[levels(SurvToth$debulking)=="suboptimal"] <- "Suboptimal"
  SurvToth$age_at_initial_pathologic_diagnosis <- SurvToth$GSE9891_eset.age_at_initial_pathologic_diagnosis
  SurvToth$vital_status <- SurvToth$GSE9891_eset.vital_status
  SurvToth$days_to_death <- SurvToth$GSE9891_eset.days_to_death
  SurvToth <- select(SurvToth, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)

#Add Yoshihara
  phenoYosh <- as.data.frame(phenoData["GSE32062.GPL6480_eset"])
  phenoYosh <- tibble::rownames_to_column(phenoYosh, "sampleid")
  SurvYosh <- merge(phenoYosh, keepall, by.x = "sampleid", by.y = "X", all.x=FALSE, all.y=FALSE)
  
  SurvYosh$stage <- SurvYosh$GSE32062.GPL6480_eset.tumorstage
  SurvYosh$debulking <- SurvYosh$GSE32062.GPL6480_eset.debulking
  levels(SurvYosh$debulking)[levels(SurvYosh$debulking)=="optimal"] <- "Optimal"
  levels(SurvYosh$debulking)[levels(SurvYosh$debulking)=="suboptimal"] <- "Suboptimal"
  SurvYosh$age_at_initial_pathologic_diagnosis <- SurvYosh$GSE32062.GPL6480_eset.age_at_initial_pathologic_diagnosis
  SurvYosh$vital_status <- SurvYosh$GSE32062.GPL6480_eset.vital_status
  SurvYosh$days_to_death <- SurvYosh$GSE32062.GPL6480_eset.days_to_death
  SurvYosh <- select(SurvYosh, sampleid, Dataset, ClusterK2_kmeans, ClusterK3_kmeans, ClusterK4_kmeans, stage, debulking, age_at_initial_pathologic_diagnosis, days_to_death, vital_status)


#Make Dataset to use for analyses

KeepCols <- c("sampleid", "Dataset", "ClusterK4_kmeans", "ClusterK3_kmeans", "age_at_initial_pathologic_diagnosis", "days_to_death", "vital_status", "debulking", "stage")

OrigTCGA <- SurvTCGA[,KeepCols]
OrigAA <- SurvAA[,KeepCols]
OrigWhite <-SurvWhite[,KeepCols]
Origmayo <- Survmayo[,KeepCols]
Origtoth <- SurvToth[,KeepCols]
Origyosh <- SurvYosh[,KeepCols]

SurvData <- rbind(OrigTCGA, OrigAA, OrigWhite, Origmayo, Origtoth, Origyosh) 


  SurvData$age[0<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<40] <- 1
  SurvData$age[40<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<45] <- 2
  SurvData$age[45<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<50] <- 3
  SurvData$age[50<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<55] <- 4
  SurvData$age[55<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<60] <- 5
  SurvData$age[60<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<65] <- 6
  SurvData$age[65<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<70] <- 7
  SurvData$age[70<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<75] <- 8
  SurvData$age[75<=SurvData$age_at_initial_pathologic_diagnosis
               & SurvData$age_at_initial_pathologic_diagnosis<110] <- 9
  
SurvData$days <- as.numeric(SurvData$days_to_death)
# Convert days into months 
SurvData$months <- as.numeric(SurvData$days) / 30.4375

SurvData$vital <- 99
SurvData$vital[SurvData$vital_status=="alive" | SurvData$vital_status=="living"] <- 0
SurvData$vital[SurvData$vital_status=="deceased"] <- 1
SurvData$vital <- as.logical(SurvData$vital)

SurvData$FewerStage <- 9
SurvData$FewerStage[SurvData$stage=="1"] <- 1
SurvData$FewerStage[SurvData$stage=="I"] <- 1
SurvData$FewerStage[SurvData$stage=="2"] <- 2
SurvData$FewerStage[SurvData$stage=="II"] <- 2
SurvData$FewerStage[SurvData$stage=="3"] <- 3
SurvData$FewerStage[SurvData$stage=="III"] <- 3
SurvData$FewerStage[SurvData$stage=="4"] <- 4
SurvData$FewerStage[SurvData$stage=="IV"] <- 4
SurvData$FewerStage[SurvData$stage=="Missing"] <- 9

levels(SurvData$debulking)[levels(SurvData$debulking)==""] <- "Missing"
levels(SurvData$debulking)[levels(SurvData$debulking)=="no surgery"] <- "Suboptimal"
levels(SurvData$debulking)[levels(SurvData$debulking)=="optimal"] <- "Optimal"
levels(SurvData$debulking)[levels(SurvData$debulking)=="suboptimal"] <- "Suboptimal"
SurvData$debulking[is.na(SurvData$debulking)]<-"Missing"

#Make Datasets to use for analyses
KeepVars <- c("sampleid", "ClusterK4_kmeans", "ClusterK3_kmeans", "age", "months", "vital", "Dataset", "debulking", "stage", "FewerStage")
KeepVarsY <- c("sampleid", "ClusterK4_kmeans", "ClusterK3_kmeans", "months", "vital", "Dataset", "debulking", "stage", "FewerStage")
AnalSet <- SurvData[,KeepVars]

#Save datasets for use in subsequent script (survival and table 1)
write.csv(AnalSet, file="AnalSet.csv")
write.csv(TableOneBlack, file="TableOneBlack.csv")
write.csv(TableOneWhite, file="TableOneWhite.csv")

