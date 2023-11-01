#This script conducts survival analyses for each of the six datasets, 
#both overall and by subtype. It also generates Kaplan-Meier curves
#to depict survival by dataset and subtype. 

library("here")
library("survival")
library("rlang")
library("dplyr")

proj_dir = here()
fig_dir = file.path(proj_dir,"figure_notebooks/survival_figs")

#readin survival datasets from survival folder data processing
AnalSetloc = file.path(proj_dir, "survival/AnalSet.csv")
AnalSet <-read.csv(AnalSetloc, header=TRUE, sep="," )
#There are 11 cases with missing stage, set them to 3, the most common stage
AnalSet$FewerStage <- ifelse (AnalSet$FewerStage==9, 3, AnalSet$FewerStage)

#Take out first 10 years of survival data
AnalSet$months <- AnalSet$months - 10
AnalSet <- AnalSet[AnalSet$months > 0, ]

#Drop observations that are completely missing
missing_id <- is.na(AnalSet$sampleid)
AnalSet <- AnalSet[!missing_id,]

SurvAACESB <-AnalSet[AnalSet$Dataset=="SchildkrautB",]
SurvAACESW <-AnalSet[AnalSet$Dataset=="SchildkrautW",]
TCGAset <-AnalSet[AnalSet$Dataset=="TCGA",]
MAYOset <-AnalSet[AnalSet$Dataset=="Mayo",]
TOTHset <-AnalSet[AnalSet$Dataset=="Tothill",]
YOSHset <-AnalSet[AnalSet$Dataset=="Yoshihara",]

SchildkrautB_n <- summary(as.factor(SurvAACESB$ClusterK4))
SchildkrautW_n <- summary(as.factor(SurvAACESW$ClusterK4))
TCGA_n <- summary(as.factor(TCGAset$ClusterK4))
MAYO_n <- summary(as.factor(MAYOset$ClusterK4))
TOTH_n <- summary(as.factor(TOTHset$ClusterK4))
YOSH_n <- summary(as.factor(YOSHset$ClusterK4))

Distributions_Drop10 <- rbind(SchildkrautB_n, SchildkrautW_n, TCGA_n,MAYO_n,TOTH_n,YOSH_n)
write.csv(Distributions_Drop10, file="distributions_drop10.csv")

###############################
##ALL
# Get a survfit object 

tmpsurvfitall <- survfit(Surv(months, vital) ~ factor(Dataset), data = AnalSet)

# Store sample size variable
nall <- nrow(AnalSet)

model <- coxph(Surv(months, vital) ~ (factor(Dataset)+ factor(age) + factor(debulking) + strata(factor(FewerStage))), AnalSet)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Wald's P
  WaldsP <- cox.sum$waldtest["pvalue"]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues, WaldsP)
  CoxPHtmploc = file.path(fig_dir, "AllData_COXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model, AnalSet)


#Likelihood ratio test of null hypothesis that there 
#is no difference in survival across datasets, adjusting for age and stage
modela <- coxph(Surv(months, vital) ~ (factor(Dataset) + factor(ClusterK4_kmeans)+ factor(FewerStage)), AnalSet)
modelb <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(FewerStage)), AnalSet)

X.lr.dataset <- 2*(modela$loglik[2]-modelb$loglik[2])
p.lr.dataset <- 1-pchisq(X.lr.dataset,5)
p.lr.dataset


##################################
##AACESB

# Get a survfit object 
tmpsurvfitAAB <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = SurvAACESB)

# Store sample size variable
nAAB <- nrow(SurvAACESB)

model3 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(age) + factor(debulking)+ strata(factor(FewerStage))), SurvAACESB)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues)
  CoxPHtmploc = file.path(fig_dir, "AACESBCOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model3, SurvAACESB)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model4 <- coxph(Surv(months, vital) ~ (factor(age)+factor(debulking)+ strata(factor(FewerStage))), SurvAACESB)

X.lr <- 2*(model3$loglik[2]-model4$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr

##################################
##AACESW

# Get a survfit object 
tmpsurvfitAAW <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = SurvAACESW)

# Store sample size variable
nAAW <- nrow(SurvAACESW)

model5 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(age) + factor(debulking)+ strata(factor(FewerStage))), SurvAACESW)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues)
  CoxPHtmploc = file.path(fig_dir, "AACESWCOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model5, SurvAACESW)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model6 <- coxph(Surv(months, vital) ~ (factor(age)+factor(debulking)+ strata(factor(FewerStage))), SurvAACESW)

X.lr <- 2*(model5$loglik[2]-model6$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr


#################################

##TCGA



# Get a survfit object 
tmpsurvfitTC <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = TCGAset)

# Store sample size variable
nTC <- nrow(TCGAset)

model7 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(age)+factor(debulking)+ strata(factor(FewerStage))), TCGAset)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues)
  CoxPHtmploc = file.path(fig_dir, "TCGACOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model7, TCGAset)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model8 <- coxph(Surv(months, vital) ~ (factor(age)+factor(debulking)+ strata(factor(FewerStage))), TCGAset)

X.lr <- 2*(model7$loglik[2]-model8$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr



#################################

##Mayo



# Get a survfit object 
tmpsurvfitMA <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = MAYOset)

# Store sample size variable
nMA <- nrow(MAYOset)

model9 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(age)+factor(debulking)+ strata(factor(FewerStage))), MAYOset)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues)
  CoxPHtmploc = file.path(fig_dir, "MAYOCOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model9, MAYOset)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model10 <- coxph(Surv(months, vital) ~ (factor(age)+factor(debulking)+ strata(factor(FewerStage))), MAYOset)

X.lr <- 2*(model9$loglik[2]-model10$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr





#################################

##Tothill



# Get a survfit object 
tmpsurvfitTO <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = TOTHset)

# Store sample size variable
nTO <- nrow(TOTHset)

model11 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(age)+factor(debulking)+ strata(factor(FewerStage))), TOTHset)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues)
  CoxPHtmploc = file.path(fig_dir, "TOTHCOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model11, TOTHset)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model12 <- coxph(Surv(months, vital) ~ (factor(age)+factor(debulking)+ strata(factor(FewerStage))), TOTHset)

X.lr <- 2*(model11$loglik[2]-model12$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr






#################################

##Yoshihara



# Get a survfit object 
tmpsurvfitYO <- survfit(Surv(months, vital) ~ ClusterK4_kmeans, data = YOSHset)

# Store sample size variable
nYO <- nrow(YOSHset)

model13 <- coxph(Surv(months, vital) ~ (factor(ClusterK4_kmeans) + factor(debulking) + strata(factor(FewerStage))), YOSHset)

coxSum <- function (coxData, name) {
  # ~~~~~~~~~~~~~~
  # This function will write out the hazards ratios, confidence intervals, pvalues,
  # and Wald's P for each cox proportional hazards model
  #
  # Args: 
  # coxData: a coxph object
  # name: the name of the dataset
  #
  # Returns:
  # a dataframe of pertinent information regarding the cox model summary
  # ~~~~~~~~~~~~~~
  
  cox.sum <- summary(coxData)
  
  # 95% Confidence Intervals
  conf.int <- cox.sum$conf.int[ ,3:4]
  
  # Hazard ratio
  hazard <- exp(coxData$coefficients)
  
  # p value
  Pvalues <- cox.sum$coefficients[ ,5]
  
  # Wald's P
  WaldsP <- cox.sum$waldtest["pvalue"]
  
  # Combine together and write to file
  CoxPHtmp <- cbind(conf.int, hazard, Pvalues, WaldsP)
  CoxPHtmploc = file.path(fig_dir, "YOSHCOXPHoutput_drop10.csv")
  write.csv(CoxPHtmp, file=CoxPHtmploc)
}

coxSum(model13, YOSHset)


#Likelihood ratio test of null hypothesis that there 
#is no association between K4 Cluster and Survival, accounting for age
model14 <- coxph(Surv(months, vital) ~ (factor(debulking)+ strata(factor(FewerStage))), YOSHset)

X.lr <- 2*(model13$loglik[2]-model14$loglik[2])
p.lr <- 1-pchisq(X.lr,3)
p.lr
