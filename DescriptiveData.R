library(table1)


CRC <- samplesheet[samplesheet$ICD10 %in% c(' C18', ' C19', ' C20'),]
Breast <- samplesheet[samplesheet$ICD10 %in% c(' C50'),]
CRC <- rbind(CRC, samplesheet[samplesheet$CanStat==F,])
Breast <- rbind(Breast, samplesheet[samplesheet$CanStat==F,])
Breast <- Breast[Breast$Sex==' F',]

Breast$BMIScore <- scale(Breast$BMIScore)
CRC$BMIScore <- scale(CRC$BMIScore)


mean(as.numeric(Breast[Breast$CanStat==T,]$Age), na.rm=T)
mean(as.numeric(Breast[Breast$CanStat==F,]$Age), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==T,]$Age), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==F,]$Age), na.rm=T)

mean(as.numeric(Breast[Breast$CanStat==T,]$Menarche), na.rm=T)
mean(as.numeric(Breast[Breast$CanStat==F,]$Menarche), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==T,]$Menarche), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==F,]$Menarche), na.rm=T)

mean(as.numeric(Breast[Breast$CanStat==T,]$BMIScore), na.rm=T)
mean(as.numeric(Breast[Breast$CanStat==F,]$BMIScore), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==T,]$BMIScore), na.rm=T)
sd(as.numeric(Breast[Breast$CanStat==F,]$BMIScore), na.rm=T)

# CRC
mean(as.numeric(CRC[CRC$CanStat==T,]$Age), na.rm=T)
mean(as.numeric(CRC[CRC$CanStat==F,]$Age), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==T,]$Age), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==F,]$Age), na.rm=T)

mean(as.numeric(CRC[CRC$CanStat==T,]$Menarche), na.rm=T)
mean(as.numeric(CRC[CRC$CanStat==F,]$Menarche), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==T,]$Menarche), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==F,]$Menarche), na.rm=T)

mean(as.numeric(CRC[CRC$CanStat==T,]$BMIScore), na.rm=T)
mean(as.numeric(CRC[CRC$CanStat==F,]$BMIScore), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==T,]$BMIScore), na.rm=T)
sd(as.numeric(CRC[CRC$CanStat==F,]$BMIScore), na.rm=T)

CRCCase <- CRC[CRC$CanStat,]
CRCControl <- CRC[-CRC$CanStat,]
sum(CRCCase$Sex==' F')/nrow(CRCCase) * 100
sum(CRCControl$Sex==' F')/nrow(CRCControl) * 100

#------------------------------------------------------------------------------#
#dbGaP data
library(readr)
library(data.table)

betas <- as.matrix(fread('C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - EpiscoresInCancer/dbGAP.txt', sep = '\t'), rownames=1)
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - EpiscoresInCancer/dbGAPEpiscores.RData")

# order beta values to be same as sample sheet
betas <- betas[,samples$SAMPLE_ID]
betatran <- t(betas)

# Set relevant variables as factors
samples$SEX <- factor(samples$SEX)
samples$RACE <- factor(samples$RACE)
samples$SMOKING <- factor(samples$SMOKING)

# Remove completely NA columns
library(tidyverse)
samples <- samples %>% select_if(~sum(!is.na(.)) > 0)
samples$BMI <- ifelse(samples$BMI<99, samples$BMI, NA)

#Marioni BMI 
Marioni_BMIMeasures <- read_csv("Marioni_BMIMeasures.csv", col_names = T, skip = 2)
Marioni_BMIMeasures <- Marioni_BMIMeasures[Marioni_BMIMeasures$CpG %in% rownames(betas),]
#Grab protein specific CpG sites
cglist <- unique(Marioni_BMIMeasures$CpG)
CGcoef <- Marioni_BMIMeasures$Coefficient

CGbeta <- betas[cglist,]

# calculate episcore
BMIScore <- CGcoef*CGbeta
if (!is.null(nrow(BMIScore))){
  BMIScore <- colSums(BMIScore)  
} 

samples$BMIScore <- BMIScore
hist(BMIScore)
rm(betas, betatran)


mean(as.numeric(samples[samples$CASE==T,]$AGE), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$AGE), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$AGE), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$AGE), na.rm=T)


mean(as.numeric(samples[samples$CASE==T,]$BMI), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$BMI), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$BMI), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$BMI), na.rm=T)

samples$BMIScore <- scale(samples$BMIScore)
mean(as.numeric(samples[samples$CASE==T,]$BMIScore), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$BMIScore), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$BMIScore), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$BMIScore), na.rm=T)



samplesCase <- samples[samples$CASE==1,]
samplesControl <- samples[samples$CASE==0,]


sum(samplesCase$SEX==1)/nrow(samplesCase) * 100
sum(samplesControl$SEX==1)/nrow(samplesControl) * 100

sum(samplesCase$RACE=='White')/nrow(samplesCase) * 100
sum(samplesControl$RACE=='White')/nrow(samplesControl) * 100


sum(samplesCase$SMOKING==3)/nrow(samplesCase) * 100
sum(samplesControl$SMOKING==3)/nrow(samplesControl) * 100

samplesCase$TTD <- samplesCase$YEAR_OF_DIAGNOSIS - samplesCase$YEAR_OF_BLOOD_DRAWN
samplesCase <- samplesCase[samplesCase$TTD <30,]
mean(samplesCase$TTD)
median(samplesCase$TTD)


#------------------------------------------------------------------------------#
# Time restricted analyses

CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- samples[samples$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

samples <- rbind(CasesOnly, ControlOnly)

mean(as.numeric(samples[samples$CASE==T,]$AGE), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$AGE), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$AGE), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$AGE), na.rm=T)


mean(as.numeric(samples[samples$CASE==T,]$BMI), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$BMI), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$BMI), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$BMI), na.rm=T)

samples$BMIScore <- scale(samples$BMIScore)
mean(as.numeric(samples[samples$CASE==T,]$BMIScore), na.rm=T)
mean(as.numeric(samples[samples$CASE==F,]$BMIScore), na.rm=T)
sd(as.numeric(samples[samples$CASE==T,]$BMIScore), na.rm=T)
sd(as.numeric(samples[samples$CASE==F,]$BMIScore), na.rm=T)





samplesCase <- samples[samples$CASE==1,]
samplesControl <- samples[samples$CASE==0,]


sum(samplesCase$SEX==1)/nrow(samplesCase) * 100
sum(samplesControl$SEX==1)/nrow(samplesControl) * 100

sum(samplesCase$RACE=='White')/nrow(samplesCase) * 100
sum(samplesControl$RACE=='White')/nrow(samplesControl) * 100


sum(samplesCase$SMOKING==3)/nrow(samplesCase) * 100
sum(samplesCase$SMOKING==2)/nrow(samplesCase) * 100
sum(samplesCase$SMOKING==1)/nrow(samplesCase) * 100
sum(samplesControl$SMOKING==3)/nrow(samplesControl) * 100

samplesCase$TTD <- samplesCase$YEAR_OF_DIAGNOSIS - samplesCase$YEAR_OF_BLOOD_DRAWN
samplesCase <- samplesCase[samplesCase$TTD <30,]
mean(samplesCase$TTD)
median(samplesCase$TTD)




