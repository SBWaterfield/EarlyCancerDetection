library(meffil)
library(GEOquery)
library(tidyverse)
library(readr)

options(mc.cores=6)

#------------------------------------------------------------------------------#
load('EPIC_Italy_GrimAge.RData')
GrimAge <- samplesheet
GrimAge$GrimAge <- resid(lm(GrimAge$GrimAge ~ GrimAge$AGE))
GrimAgeResid <- GrimAge$GrimAge
  

load('EPIC_Italy_Meffil.RData')
load('EPIC_Italy_Samples.RData')

beta <- beta$beta
samplesheet <- samplesheet[samplesheet$Sample_Name %in% colnames(beta),]
betas <-t(beta)



# Build PC clocks
samplesheet$Female <- ifelse(samplesheet$Sex==' F', 1, 0)
samplesheet$Age <- as.numeric(samplesheet$Age)

# clocksDir <- getwd()
# source(paste(clocksDir, "/run_calcPCClocks.R", sep = ""))
# source(paste(clocksDir, "/run_calcPCClocks_Accel.R", sep = ""))
# clocksDir <- paste(clocksDir, "/", sep = "")
# PCClock_DNAmAge <- calcPCClocks(path_to_PCClocks_directory = clocksDir, datMeth = betatran, datPheno = samplesheet)
# PCClock_DNAmAgeAccel <- calcPCClocks_Accel(PCClock_DNAmAge)
# save(PCClock_DNAmAge, PCClock_DNAmAgeAccel, file = 'PCClock_EPIC_Italy.RData')

load('PCClock_EPIC_Italy.RData')
samplesheet <- cbind(samplesheet, PCClock_DNAmAgeAccel)

load('EPICItaly_StandardClocks.RData')
colnames(acc) <- paste0(colnames(acc), 'Resid')
samplesheet <- cbind(samplesheet, acc)

samplesheet <- cbind(samplesheet, GrimAgeResid)


# Load in Gadd Coefficients
Episcore_Coef <- read_csv("Data/Episcore_Coef.csv")



# Build the episcores
Soma <-  Episcore_Coef[Episcore_Coef$Panel == 'SomaScan',]
Olink <- Episcore_Coef[Episcore_Coef$Panel == 'Olink',]
#ArrayIntersect <- intersect(Soma$`Gene Name`, Olink$`Gene Name`) # 5 shared proteins

# Specify Olink proteins
Olink$`Gene Name` <- paste(Olink$`Gene Name`, '_Olink', sep = '')

# Combine dataframes
Episcore_Coef <- rbind(Soma, Olink)
Episcore_Genes <- unique(Episcore_Coef$`Gene Name`) #108, 5 of which are in both soma/olink

# Reduce Episcore_Coefs to CpGs in beta values
betanames <- rownames(beta)
Episcore_Coef <- Episcore_Coef[Episcore_Coef$`CpG Site` %in% betanames,] 

#------------------------------------------------------------------------------#
# Empty table
Episcores <- data.frame()

# Loop to create episcores in siNET data
for(Protein in Episcore_Genes){
  
  # Load in Protein specific episcore details
  df <- data.frame(Episcore_Coef[Episcore_Coef$`Gene Name`==Protein,])  
  
  #Grab protein specific CpG sites
  cglist <- unique(df$CpG.Site)
  CGcoef <- df$CpG.Coeficient
  CGbeta <- beta[cglist,]
  
  # calculate episcore
  Episcore <- CGcoef*CGbeta
  if (!is.null(nrow(Episcore))){
    Episcore <- colSums(Episcore)
  }
  # Add to Table of data
  newrow <- c(Protein, Episcore)
  Episcores <- rbind(Episcores, newrow)
}

# Column names
Episcore_Cols <- c('Protein')
Episcore_Cols <- append(Episcore_Cols, colnames(beta))
colnames(Episcores) <- Episcore_Cols

# Make episcores gene name column the row name
Episcores <- data.frame(Episcores[,-1], row.names = Episcores[,1])
i <- c(1:ncol(beta))    
Episcores[ , i] <- apply(Episcores[ , i], 2, function(x) as.numeric(as.character(x)))
Episcores <- t(Episcores)

# Bind Episcores to samplesheetheet 
samplesheet <- cbind(samplesheet, Episcores)
#------------------------------------------------------------------------------#
#Marioni BMI 
Marioni_BMIMeasures <- read_csv("Marioni_BMIMeasures.csv", col_names = T, skip = 2)

#Grab protein specific CpG sites
cglist <- unique(Marioni_BMIMeasures$CpG)
CGcoef <- Marioni_BMIMeasures$Coefficient
CGbeta <- beta[cglist,]

# calculate episcore
BMIScore <- CGcoef*CGbeta
if (!is.null(nrow(BMIScore))){
  BMIScore <- colSums(BMIScore)  
}  
samplesheet <- cbind(samplesheet, BMIScore)

samplesheetplaceholder <- samplesheet

#save(samplesheetplaceholder ,file = 'EPICItalysamplesheetplaceholder.RData')
#------------------------------------------------------------------------------#

EpiscoreNames <- Episcore_Genes
ClockNames <- c('PCHorvath1Resid', 'PCHorvath2Resid', 'PCHannumResid', 'PCPhenoAgeResid', 'PCDNAmTLResid', 'PCGrimAgeResid',
                'horvathResid', 'hannumResid', 'phenoageResid', 'dunedinpaceResid', 'dunedinpoam38Resid', 'brennerResid', 'dnamtlResid', 'GrimAgeResid')


#------------------------------------------------------------------------------#
samplesheet <- samplesheetplaceholder
samplesheet$CanStat <- as.factor(samplesheet$CanStat)
unique(samplesheet$ICD10)
samplesheet$ICD10[is.na(samplesheet$ICD10)] <- 'Control' 
samplesheet$Age <- as.numeric(samplesheet$Age)
samplesheet$Menarche <- as.numeric(samplesheet$Menarche)
samplesheet$TTD <- as.numeric(samplesheet$TTD)

samplesheet <- samplesheet[samplesheet$ICD10 == ' C50',]
samplesheet <- samplesheet[samplesheet$Sex ==' F',]
EPICBreastTTDReg <- data.frame()
for (Episcore in EpiscoreNames) {
  
  samplesheet[[Episcore]] <- scale(samplesheet[[Episcore]])
  
  #model <- glm(CanStat ~ samplesheet[[Episcore]], family = 'binomial', data=samplesheet)
  model <- glm(TTD ~ samplesheet[[Episcore]] + Age + Menarche +BMIScore, family = 'gaussian', data=samplesheet)
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(Episcore, modsum)
  EPICBreastTTDReg <- rbind(EPICBreastTTDReg, newrow)
  
  
  
  
  
}
colnames(EPICBreastTTDReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')

EPICBreastTTDReg$Coef <- as.numeric(EPICBreastTTDReg$Coef)
EPICBreastTTDReg$SE <- as.numeric(EPICBreastTTDReg$SE)

EPICBreastTTDReg$OR <- exp(EPICBreastTTDReg$Coef)
EPICBreastTTDReg$LI <- exp(EPICBreastTTDReg$Coef-1.96*EPICBreastTTDReg$SE)
EPICBreastTTDReg$UI <- exp(EPICBreastTTDReg$Coef+1.96*EPICBreastTTDReg$SE)


EPICBreastTTDReg_Gadd <- EPICBreastTTDReg[!str_detect(EPICBreastTTDReg$Episcore, 'Resid'),]
EPICBreastTTDReg_Gadd$FDR <- p.adjust(EPICBreastTTDReg_Gadd$P, method = 'fdr')
write.csv(EPICBreastTTDReg_Gadd, file = 'EPICBreast_EpiscoreGadd_TTD.csv')

EPICBreastTTDReg <- data.frame()
for (Episcore in ClockNames) {
  
  samplesheet[[Episcore]] <- scale(samplesheet[[Episcore]])
  
  #model <- glm(CanStat ~ samplesheet[[Episcore]], family = 'binomial', data=samplesheet)
  model <- glm(TTD ~ samplesheet[[Episcore]] + Menarche + BMIScore, family = 'gaussian', data=samplesheet)
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(Episcore, modsum)
  EPICBreastTTDReg <- rbind(EPICBreastTTDReg, newrow)
  
  
  
  
  
}
colnames(EPICBreastTTDReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')

EPICBreastTTDReg$Coef <- as.numeric(EPICBreastTTDReg$Coef)
EPICBreastTTDReg$SE <- as.numeric(EPICBreastTTDReg$SE)
EPICBreastTTDReg$OR <- exp(EPICBreastTTDReg$Coef)
EPICBreastTTDReg$LI <- exp(EPICBreastTTDReg$Coef-1.96*EPICBreastTTDReg$SE)
EPICBreastTTDReg$UI <- exp(EPICBreastTTDReg$Coef+1.96*EPICBreastTTDReg$SE)
EPICBreastTTDReg_Clock <- EPICBreastTTDReg[str_detect(EPICBreastTTDReg$Episcore, 'Resid'),]
EPICBreastTTDReg_Clock$FDR <- p.adjust(EPICBreastTTDReg_Clock$P, method = 'fdr')
write.csv(EPICBreastTTDReg_Clock, file = 'EPICBreast_EpiscoreClock_TTD.csv')

#------------------------------------------------------------------------------#
# Cox PH models episcores
library(survival)
library(survminer)

samples <- samplesheetplaceholder
samples$ICD10[is.na(samples$ICD10)] <- 'Control'
samples$Sex <- as.factor(samples$Sex)
samples$Age <- as.numeric(samples$Age)
samples <- samples[samples$Sex ==' F',]
samples$Menarche <- as.numeric(samples$Menarche)
samples$TTD <- as.numeric(samples$TTD)
CasesOnly <- samples[samples$ICD10==' C50',]
ControlOnly <- samples[samples$ICD10=='Control',]
ControlOnly$TTD <- max(samples$TTD, na.rm = T)

coxdata <- rbind(CasesOnly, ControlOnly)
coxdata$CanStat <- ifelse(coxdata$CanStat==T, 1,0)

EpiscoreCox_TTDiag <- data.frame()

for (episcore in EpiscoreNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  res.cox <- coxph(Surv(TTD, CanStat) ~ episcore +Age + Menarche + BMIScore, data = coxdata)
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ IGFBP1 + SEX + AGE + BMI + SMOKING, data = coxdata)
  
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  

}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag[!str_detect(EpiscoreCox_TTDiag$Episcore, 'Resid'),]


EpiscoreCox_TTDiag_Gadd$FDR <- p.adjust(EpiscoreCox_TTDiag_Gadd$P, method = 'fdr')


EpiscoreCox_TTDiag_Gadd[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Gadd[[x]]))
EpiscoreCox_TTDiag_Gadd[,2:7] <- signif(EpiscoreCox_TTDiag_Gadd[,2:7],3)
write.csv(EpiscoreCox_TTDiag_Gadd, file = 'EPICBreast_EpiscoreGadd_Cox.csv')


EpiscoreCox_TTDiag <- data.frame()

for (episcore in ClockNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  res.cox <- coxph(Surv(TTD, CanStat) ~ episcore + Menarche + BMIScore, data = coxdata)
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ IGFBP1 + SEX + AGE + BMI + SMOKING, data = coxdata)
  
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  

  
}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')
EpiscoreCox_TTDiag_Clock <- EpiscoreCox_TTDiag[str_detect(EpiscoreCox_TTDiag$Episcore, 'Resid'),]


EpiscoreCox_TTDiag_Clock$FDR <- p.adjust(EpiscoreCox_TTDiag_Clock$P, method = 'fdr')


EpiscoreCox_TTDiag_Clock[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Clock[[x]]))
EpiscoreCox_TTDiag_Clock[,2:7] <- signif(EpiscoreCox_TTDiag_Clock[,2:7],3)
write.csv(EpiscoreCox_TTDiag_Clock, file = 'EPICBreast_EpiscoreClock_Cox.csv')
#------------------------------------------------------------------------------#
# OSCA analyses
samplesheet <- samplesheetplaceholder
samplesheet$CanStat <- as.factor(samplesheet$CanStat)
unique(samplesheet$ICD10)
samplesheet$ICD10[is.na(samplesheet$ICD10)] <- 'Control' 

samplesheet$Sex <- ifelse(samplesheet$Sex == ' F', 'F', 'M')
samplesheet$Sex <- as.factor(samplesheet$Sex)
samplesheet$Age <- as.numeric(samplesheet$Age)
samplesheet$Menarche <- as.numeric(samplesheet$Menarche)
samplesheet$TTD <- as.numeric(samplesheet$TTD)

samplesheet <- samplesheet[samplesheet$ICD10 %in% c(' C50'),]

for (episcore in Episcore_Genes){
  
  res <- resid(lm(samplesheet[[episcore]] ~ samplesheet$Age))
  samplesheet[[episcore]] <- res
}


OSCAFun <- function(df, XVar, Var, Resname = 'EpiscoresBreast'){
  
  df$IID <- df$Sample_Name
  df$FID <- df$Sample_Name
  dfplaceholder <- df
  #Make ORM files
  df <- df
  
  df <- df[c('FID', 'IID', XVar)]
  
  df <- df[!is.na(df$IID),]
  
  df <- df %>% mutate_at(c(XVar), ~(scale(.) %>% as.vector))
  
  
  write.table(df, file = 'temptxt.txt', sep=' ', row.names = F, col.names = T)
  
  system('./osca-0.46.1 --efile temptxt.txt --make-bod --out myprofile')
  system('./osca-0.46.1 --befile myprofile --make-orm --out myorm')
  
  # name <- Var
  # file.rename(from='myorm.orm.N.bin', to= paste(name, '.orm.N.bin', sep = ''))
  # file.rename(from='myorm.orm.id', to= paste(name, '.orm.id', sep = ''))
  # file.rename(from='myorm.orm.bin', to= paste(name, '.orm.bin', sep = ''))
  # 
  # #Change myprofile name for matching  of data in MOA stage
  # file.rename(from='myprofile.bod', to= paste(name, '.bod', sep = ''))
  # file.rename(from='myprofile.opi', to= paste(name, '.opi', sep = ''))
  # file.rename(from='myprofile.oii', to= paste(name, '.oii', sep = ''))
  
  df <- dfplaceholder
  SexCovar <- df[c('FID', 'IID', 'Sex')]
  CellCovar <- df[c('FID', 'IID', 'BMIScore', 'Menarche')]
  
  
  df <- df[c('FID', 'IID', Var)]
  df <- df[!is.na(df$IID),]
  colnames(df) <- c('FID', 'IID', 'pheno')
  
  
  
  write.table(df, file ='phenotest.txt',  sep=' ', row.names = F, col.names = T)
  write.table(SexCovar, file ='SexCovar.txt',  sep=' ', row.names = F, col.names = T)
  write.table(CellCovar, file ='CellCovar.txt',  sep=' ', row.names = F, col.names = T)
  
  
  #Run REML
  system('./osca-0.46.1 --reml --orm myorm --pheno phenotest.txt --qcovar CellCovar.txt --out myreml')
  
  
  # Run MOA 
  #system(paste0(
  #  './osca-0.46.1 --moa-exact --befile ', name, ' --pheno phenotest.txt --covar SexCovar.txt  --qcovar CellCovar.txt --orm ', name, ' --out myMOA'
  #))
  
  
  # Manage the weird data output
  read_the_text <- scan(file = "myreml.rsq", what = character(),sep = "\n")
  
  split_each_line_by_spaces <- strsplit(x = read_the_text,split = "\t")
  
  get_element_names <- lapply(X = split_each_line_by_spaces,FUN = `[[`,i = 1)
  
  get_element_values <- lapply(X = split_each_line_by_spaces,FUN = `[`,i = (-1))
  
  required_result_as_character <- setNames(object = get_element_values,nm = get_element_names)
  
  required_result <- lapply(X = required_result_as_character,FUN = as.numeric)
  
  # Grab results of REML
  newresults <- c(required_result[["V(O)/Vp"]][1], required_result[["Pval"]], 
                  required_result[["n"]], Resname)
  
  file.rename(from='myreml.rsq', to= paste(Resname, '_', 'reml.rsq', sep = ''))
  
  return(newresults)

} 

ClocksNames <-  setdiff(EpiscoreNames, Episcore_Genes)
# EpiscoresREML <- OSCAFun(samplesheet, Episcore_Genes, 'TTD', 'Episcores_Breast')
# EpiscoresClocksREML <- OSCAFun(samplesheet, EpiscoreNames, 'TTD', 'Episcores_and_Clocks_Breast')
# ClocksREML <- OSCAFun(samplesheet, ClocksNames, 'TTD', 'Clocks_Breast')

#OSCA_Breast_Res <- rbind(EpiscoresREML, EpiscoresClocksREML, ClocksREML)
#colnames(OSCA_Breast_Res) <- c('V(0)/Vp', 'P', 'SampleN', 'Analysis')
#save(OSCA_Breast_Res, file = 'OSCA_Breast_TTD.RData')

OSCAFun(samplesheet, Episcore_Genes, 'TTD', 'Episcores_Breast')
OSCAFun(samplesheet, EpiscoreNames, 'TTD', 'Episcores_and_Clocks_Breast')
OSCAFun(samplesheet, ClocksNames, 'TTD', 'Clocks_Breast')

OSCAFun <- function(df, XVar, Var, Resname = 'Breast'){
  
  df$IID <- df$Sample_Name
  df$FID <- df$Sample_Name
  dfplaceholder <- df
  #Make ORM files
  df <- df
  
  df <- df[c('FID', 'IID')]
  
  samplenames <- samplesheet$Sample_Name
  XVar <- XVar[samplenames,]

  df <- cbind(df, XVar)
  
  df <- df[!is.na(df$IID),]
  
  #df <- df %>% mutate_at(c(XVar), ~(scale(.) %>% as.vector))
  
  
  write.table(df, file = 'temptxt.txt', sep=' ', row.names = F, col.names = T)
  
  system('./osca-0.46.1 --efile temptxt.txt --make-bod --out myprofile')
  system('./osca-0.46.1 --befile myprofile --make-orm --out myorm')
  
  # name <- Var
  # file.rename(from='myorm.orm.N.bin', to= paste(name, '.orm.N.bin', sep = ''))
  # file.rename(from='myorm.orm.id', to= paste(name, '.orm.id', sep = ''))
  # file.rename(from='myorm.orm.bin', to= paste(name, '.orm.bin', sep = ''))
  # 
  # #Change myprofile name for matching  of data in MOA stage
  # file.rename(from='myprofile.bod', to= paste(name, '.bod', sep = ''))
  # file.rename(from='myprofile.opi', to= paste(name, '.opi', sep = ''))
  # file.rename(from='myprofile.oii', to= paste(name, '.oii', sep = ''))
  
  df <- dfplaceholder
  SexCovar <- df[c('FID', 'IID', 'Sex')]
  CellCovar <- df[c('FID', 'IID', 'BMIScore', 'Menarche')]
  
  
  df <- df[c('FID', 'IID', Var)]
  df <- df[!is.na(df$IID),]
  colnames(df) <- c('FID', 'IID', 'pheno')
  
  write.table(df, file ='phenotest.txt',  sep=' ', row.names = F, col.names = T)
  write.table(SexCovar, file ='SexCovar.txt',  sep=' ', row.names = F, col.names = T)
  write.table(CellCovar, file ='CellCovar.txt',  sep=' ', row.names = F, col.names = T)
  
  
  #Run REML
  system('./osca-0.46.1 --reml --orm myorm --pheno phenotest.txt --qcovar CellCovar.txt --out myreml')
  
  
  # Run MOA 
  #system(paste0(
  #  './osca-0.46.1 --moa-exact --befile ', name, ' --pheno phenotest.txt --covar SexCovar.txt  --qcovar CellCovar.txt --orm ', name, ' --out myMOA'
  #))
  
  
  # Manage the weird data output
  read_the_text <- scan(file = "myreml.rsq", what = character(),sep = "\n")
  
  split_each_line_by_spaces <- strsplit(x = read_the_text,split = "\t")
  
  get_element_names <- lapply(X = split_each_line_by_spaces,FUN = `[[`,i = 1)
  
  get_element_values <- lapply(X = split_each_line_by_spaces,FUN = `[`,i = (-1))
  
  required_result_as_character <- setNames(object = get_element_values,nm = get_element_names)
  
  required_result <- lapply(X = required_result_as_character,FUN = as.numeric)
  
  # Grab results of REML
  newresults <- c(required_result[["V(O)/Vp"]][1], required_result[["Pval"]], 
                  required_result[["n"]], Resname)
  
  file.rename(from='myreml.rsq', to= paste(Resname, '_', 'reml.rsq', sep = ''))
  
  return(newresults)

} 


#OSCAEWAS_Breast_Res <- OSCAFun(samplesheet, betas, 'TTD', 'EWAS_Breast')
#colnames(OSCAEWAS_Breast_Res) <- c('V(0)/Vp', 'P', 'SampleN', 'Analysis')
#save(OSCAEWAS_Breast_Res, file = 'OSCAEWAS_Breast_TTD.RData')

OSCAFun(samplesheet, betas, 'TTD', 'EWAS_Breast')
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#-------------------------Colorectal Cancer------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#TTD in colrectal cancer 
samplesheet <- samplesheetplaceholder
samplesheet$CanStat <- as.factor(samplesheet$CanStat)
unique(samplesheet$ICD10)
samplesheet$ICD10[is.na(samplesheet$ICD10)] <- 'Control' 
samplesheet$Sex <- as.factor(samplesheet$Sex)
samplesheet$Age <- as.numeric(samplesheet$Age)
samplesheet$Menarche <- as.numeric(samplesheet$Menarche)
samplesheet$TTD <- as.numeric(samplesheet$TTD)

samplesheet <- samplesheet[samplesheet$ICD10 %in% c(' C18', ' C20', ' C19'),]

EPICCCRTTDReg <- data.frame()
for (Episcore in EpiscoreNames) {
  
  samplesheet[[Episcore]] <- scale(samplesheet[[Episcore]])
  #model <- glm(CanStat ~ samplesheet[[Episcore]], family = 'binomial', data=samplesheet)
  model <- glm(TTD ~ samplesheet[[Episcore]] + Age +Sex, family = 'gaussian', data=samplesheet)
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(Episcore, modsum)
  EPICCCRTTDReg <- rbind(EPICCCRTTDReg, newrow)
  
  
  
  
  
}
colnames(EPICCCRTTDReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')

EPICCCRTTDReg_Gadd <- EPICCCRTTDReg[!str_detect(EPICCCRTTDReg$Episcore, 'Resid'),]
EPICCCRTTDReg_Gadd$FDR <- p.adjust(EPICCCRTTDReg_Gadd$P, method = 'fdr')
write.csv(EPICCCRTTDReg_Gadd, file = 'EPICCCR_EpiscoreGadd_TTD.csv')

EPICCCRTTDReg <- data.frame()
for (Episcore in ClockNames) {
  
  samplesheet[[Episcore]] <- scale(samplesheet[[Episcore]])
  #model <- glm(CanStat ~ samplesheet[[Episcore]], family = 'binomial', data=samplesheet)
  model <- glm(TTD ~ samplesheet[[Episcore]] + Sex, family = 'gaussian', data=samplesheet)
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(Episcore, modsum)
  EPICCCRTTDReg <- rbind(EPICCCRTTDReg, newrow)
  
  
  
  
  
}
colnames(EPICCCRTTDReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')

EPICCCRTTDReg_Clock <- EPICCCRTTDReg[str_detect(EPICCCRTTDReg$Episcore, 'Resid'),]
EPICCCRTTDReg_Clock$FDR <- p.adjust(EPICCCRTTDReg_Clock$P, method = 'fdr')
write.csv(EPICCCRTTDReg_Clock, file = 'EPICCCR_EpiscoreClock_TTD.csv')
#==============================================================================#
# Cox PH models 
library(survival)
library(survminer)

samples <- samplesheetplaceholder
samples$ICD10[is.na(samples$ICD10)] <- 'Control'
samples$Sex <- as.factor(samples$Sex)
samples$Age <- as.numeric(samples$Age)
samples$Menarche <- as.numeric(samples$Menarche)
samples$TTD <- as.numeric(samples$TTD)
CasesOnly <- samples[samples$ICD10 %in% c(' C18', ' C19', ' C20'),]

ControlOnly <- samples[samples$ICD10=='Control',]
ControlOnly$TTD <- max(CasesOnly$TTD)

coxdata <- rbind(CasesOnly, ControlOnly)
coxdata$CanStat <- ifelse(coxdata$CanStat==T, 1,0)

EpiscoreCox_TTDiag <- data.frame()

for (episcore in EpiscoreNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  
  res.cox <- coxph(Surv(TTD, CanStat) ~ episcore + Age + Sex + BMIScore, data = coxdata)
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ IGFBP1 + SEX + AGE + BMI + SMOKING, data = coxdata)
  
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  

  
}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')


EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag[!str_detect(EpiscoreCox_TTDiag$Episcore, 'Resid'),]
EpiscoreCox_TTDiag_Gadd$FDR <- p.adjust(EpiscoreCox_TTDiag_Gadd$P, method = 'fdr')
EpiscoreCox_TTDiag_Gadd[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Gadd[[x]]))
EpiscoreCox_TTDiag_Gadd[,2:7] <- signif(EpiscoreCox_TTDiag_Gadd[,2:7],3)
write.csv(EpiscoreCox_TTDiag_Gadd, file = 'EPICCCR_EpiscoreGadd_Cox.csv')


EpiscoreCox_TTDiag <- data.frame()

for (episcore in ClockNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  
  res.cox <- coxph(Surv(TTD, CanStat) ~ episcore + Sex + BMIScore, data = coxdata)
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ IGFBP1 + SEX + AGE + BMI + SMOKING, data = coxdata)
  
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  

}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')


EpiscoreCox_TTDiag_Clock <- EpiscoreCox_TTDiag[str_detect(EpiscoreCox_TTDiag$Episcore, 'Resid'),]
EpiscoreCox_TTDiag_Clock$FDR <- p.adjust(EpiscoreCox_TTDiag_Clock$P, method = 'fdr')
EpiscoreCox_TTDiag_Clock[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Clock[[x]]))
EpiscoreCox_TTDiag_Clock[,2:7] <- signif(EpiscoreCox_TTDiag_Clock[,2:7],3)
write.csv(EpiscoreCox_TTDiag_Clock, file = 'EPICCCR_EpiscoreClock_Cox.csv')
#------------------------------------------------------------------------------#
# OSCA analyses
samplesheet <- samplesheetplaceholder
samplesheet$CanStat <- as.factor(samplesheet$CanStat)
unique(samplesheet$ICD10)
samplesheet$ICD10[is.na(samplesheet$ICD10)] <- 'Control' 

samplesheet$Sex <- ifelse(samplesheet$Sex == ' F', 'F', 'M')
samplesheet$Sex <- as.factor(samplesheet$Sex)
samplesheet$Age <- as.numeric(samplesheet$Age)
samplesheet$Menarche <- as.numeric(samplesheet$Menarche)
samplesheet$TTD <- as.numeric(samplesheet$TTD)

samplesheet <- samplesheet[samplesheet$ICD10 %in% c(' C18', ' C20', ' C19'),]


for (episcore in Episcore_Genes){
  
  res <- resid(lm(samplesheet[[episcore]] ~ samplesheet$Age))
  samplesheet[[episcore]] <- res
}


  

OSCAFun <- function(df, XVar, Var, Resname = 'EpiscoresCRC'){
  
  df$IID <- df$Sample_Name
  df$FID <- df$Sample_Name
  dfplaceholder <- df
  #Make ORM files
  df <- df

  df <- df[c('FID', 'IID', XVar)]
  
  df <- df[!is.na(df$IID),]
  
  df <- df %>% mutate_at(c(XVar), ~(scale(.) %>% as.vector))
  
  
  write.table(df, file = 'temptxt.txt', sep=' ', row.names = F, col.names = T)
  
  system('./osca-0.46.1 --efile temptxt.txt --make-bod --out myprofile')
  system('./osca-0.46.1 --befile myprofile --make-orm --out myorm')

  # name <- Var
  # file.rename(from='myorm.orm.N.bin', to= paste(name, '.orm.N.bin', sep = ''))
  # file.rename(from='myorm.orm.id', to= paste(name, '.orm.id', sep = ''))
  # file.rename(from='myorm.orm.bin', to= paste(name, '.orm.bin', sep = ''))
  # 
  # #Change myprofile name for matching  of data in MOA stage
  # file.rename(from='myprofile.bod', to= paste(name, '.bod', sep = ''))
  # file.rename(from='myprofile.opi', to= paste(name, '.opi', sep = ''))
  # file.rename(from='myprofile.oii', to= paste(name, '.oii', sep = ''))
  
  df <- dfplaceholder
  SexCovar <- df[c('FID', 'IID', 'Sex')]
  CellCovar <- df[c('FID', 'IID', 'BMIScore')]

  
  df <- df[c('FID', 'IID', Var)]
  df <- df[!is.na(df$IID),]
  colnames(df) <- c('FID', 'IID', 'pheno')
  
  
  
  write.table(df, file ='phenotest.txt',  sep=' ', row.names = F, col.names = T)
  write.table(SexCovar, file ='SexCovar.txt',  sep=' ', row.names = F, col.names = T)
  write.table(CellCovar, file ='CellCovar.txt',  sep=' ', row.names = F, col.names = T)
  
  
  #Run REML
  system('./osca-0.46.1 --reml --orm myorm --pheno phenotest.txt --covar SexCovar.txt --qcovar CellCovar.txt --out myreml')

  
  # Run MOA 
  #system(paste0(
  #  './osca-0.46.1 --moa-exact --befile ', name, ' --pheno phenotest.txt --covar SexCovar.txt  --qcovar CellCovar.txt --orm ', name, ' --out myMOA'
  #))
  
  
  # Manage the weird data output
  read_the_text <- scan(file = "myreml.rsq", what = character(),sep = "\n")
  
  split_each_line_by_spaces <- strsplit(x = read_the_text,split = "\t")
  
  get_element_names <- lapply(X = split_each_line_by_spaces,FUN = `[[`,i = 1)
  
  get_element_values <- lapply(X = split_each_line_by_spaces,FUN = `[`,i = (-1))
  
  required_result_as_character <- setNames(object = get_element_values,nm = get_element_names)
  
  required_result <- lapply(X = required_result_as_character,FUN = as.numeric)
  
  # Grab results of REML
  newresults <- c(required_result[["V(O)/Vp"]][1], required_result[["Pval"]], 
                  required_result[["n"]], Resname)

  file.rename(from='myreml.rsq', to= paste(Resname, '_', 'reml.rsq', sep = ''))
  
  return(newresults)

} 

ClocksNames <-  setdiff(EpiscoreNames, Episcore_Genes)
# EpiscoresREML <- OSCAFun(samplesheet, Episcore_Genes, 'TTD', 'Episcores_CRC')
# EpiscoresClocksREML <- OSCAFun(samplesheet, EpiscoreNames, 'TTD', 'Episcores_and_Clocks_CRC')
# ClocksREML <- OSCAFun(samplesheet, ClocksNames, 'TTD', 'Clocks_CRC')
# 
# OSCA_CRC_Res <- rbind(EpiscoresREML, EpiscoresClocksREML, ClocksREML)
# colnames(OSCA_CRC_Res) <- c('V(0)/Vp', 'P', 'SampleN', 'Analysis')
# save(OSCA_CRC_Res, file = 'OSCA_CRC_TTD.RData')

OSCAFun(samplesheet, Episcore_Genes, 'TTD', 'Episcores_CRC')
OSCAFun(samplesheet, EpiscoreNames, 'TTD', 'Episcores_and_Clocks_CRC')
OSCAFun(samplesheet, ClocksNames, 'TTD', 'Clocks_CRC')
#------------------------------------------------------------------------------#

# Methylation REML


OSCAFun <- function(df, XVar, Var, Resname = 'CRC'){
  
  df$IID <- df$Sample_Name
  df$FID <- df$Sample_Name
  dfplaceholder <- df
  #Make ORM files
  df <- df
  
  df <- df[c('FID', 'IID')]
  
  samplenames <- samplesheet$Sample_Name
  XVar <- XVar[samplenames,]
  
  df <- cbind(df, XVar)
  
  df <- df[!is.na(df$IID),]
  
  #df <- df %>% mutate_at(c(XVar), ~(scale(.) %>% as.vector))
  
  
  write.table(df, file = 'temptxt.txt', sep=' ', row.names = F, col.names = T)
  
  system('./osca-0.46.1 --efile temptxt.txt --make-bod --out myprofile')
  system('./osca-0.46.1 --befile myprofile --make-orm --out myorm')
  
  # name <- Var
  # file.rename(from='myorm.orm.N.bin', to= paste(name, '.orm.N.bin', sep = ''))
  # file.rename(from='myorm.orm.id', to= paste(name, '.orm.id', sep = ''))
  # file.rename(from='myorm.orm.bin', to= paste(name, '.orm.bin', sep = ''))
  # 
  # #Change myprofile name for matching  of data in MOA stage
  # file.rename(from='myprofile.bod', to= paste(name, '.bod', sep = ''))
  # file.rename(from='myprofile.opi', to= paste(name, '.opi', sep = ''))
  # file.rename(from='myprofile.oii', to= paste(name, '.oii', sep = ''))
  
  df <- dfplaceholder
  SexCovar <- df[c('FID', 'IID', 'Sex')]
  CellCovar <- df[c('FID', 'IID', 'BMIScore')]
  
  
  df <- df[c('FID', 'IID', Var)]
  df <- df[!is.na(df$IID),]
  colnames(df) <- c('FID', 'IID', 'pheno')
  
  write.table(df, file ='phenotest.txt',  sep=' ', row.names = F, col.names = T)
  write.table(SexCovar, file ='SexCovar.txt',  sep=' ', row.names = F, col.names = T)
  write.table(CellCovar, file ='CellCovar.txt',  sep=' ', row.names = F, col.names = T)
  
  
  #Run REML
  system('./osca-0.46.1 --reml --orm myorm --pheno phenotest.txt --covar SexCovar.txt --qcovar CellCovar.txt --out myreml')
  
  
  # Run MOA 
  #system(paste0(
  #  './osca-0.46.1 --moa-exact --befile ', name, ' --pheno phenotest.txt --covar SexCovar.txt  --qcovar CellCovar.txt --orm ', name, ' --out myMOA'
  #))
  
  
  # Manage the weird data output
  read_the_text <- scan(file = "myreml.rsq", what = character(),sep = "\n")
  
  split_each_line_by_spaces <- strsplit(x = read_the_text,split = "\t")
  
  get_element_names <- lapply(X = split_each_line_by_spaces,FUN = `[[`,i = 1)
  
  get_element_values <- lapply(X = split_each_line_by_spaces,FUN = `[`,i = (-1))
  
  required_result_as_character <- setNames(object = get_element_values,nm = get_element_names)
  
  required_result <- lapply(X = required_result_as_character,FUN = as.numeric)
  
  # Grab results of REML
  newresults <- c(required_result[["V(O)/Vp"]][1], required_result[["Pval"]], 
                  required_result[["n"]], Resname)
  
  file.rename(from='myreml.rsq', to= paste(Resname, '_', 'reml.rsq', sep = ''))
  return(newresults)
  

} 


#OSCAEWAS_CRC_Res <- OSCAFun(samplesheet, betas, 'TTD', 'EWAS_CRC')
#save(OSCAEWAS_Panc_Res, file = 'OSCAEWAS_CRC_TTD.RData')

OSCAFun(samplesheet, betas, 'TTD', 'EWAS_CRC')





