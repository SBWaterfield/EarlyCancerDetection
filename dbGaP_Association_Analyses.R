library(InformationValue)
library(caret)
library(ggplot2)
library(lme4)
library(dplyr)
library(pROC)
library(stringr)
library(randomForest)
library(glmnet)
library(glmnetUtils)

load("dbGAPEpiscores.RData")


# Set relevant variables as factors
samples$SEX <- factor(samples$SEX)
samples$RACE <- factor(samples$RACE)
samples$SMOKING <- factor(samples$SMOKING)

# Remove completely NA columns
samples <- samples %>% select_if(~sum(!is.na(.)) > 0)
samples$BMI <- ifelse(samples$BMI<99, samples$BMI, NA)

samplesplaceholder <- samples


EpiscoreNames <- colnames(samples)[17:ncol(samples)]

load('dbGAPCellCounts.RData')
samplesplaceholder <- cbind(samplesplaceholder, CellCounts)

load("dbGAPClocks.RData")
samples<- samples[,c('horvath', 'hannum', 'phenoage', 'brenner', 'dnamtl', 'GrimAge')]

load('dbclocks.RData')
dbclocks <- dbclocks[, c('dunedinpace', 'dunedinpoam38')]

samplesplaceholder <- cbind(samplesplaceholder, samples, dbclocks)

#------------------------------------------------------------------------------#
#-------------------------Time to Diagnosis------------------------------------#
#------------------------------------------------------------------------------#
# Test basic association between time to diagnosis and episcores
samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=30,]

EpiscoreDiagnosticGap <- data.frame()

for (episcore in EpiscoreNames) {
  
  CasesOnly[[episcore]] <- scale(CasesOnly[[episcore]])
  model <- glm(TimeGap ~ CasesOnly[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(TimeGap ~ CasesOnly[[episcore]], family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreDiagnosticGap <- rbind(EpiscoreDiagnosticGap, newrow)
  
}
colnames(EpiscoreDiagnosticGap) <- c('Episcore', 'Coef', 'SE', 'B', 'P')

EpiscoreDiagnosticGap_Gadd <- EpiscoreDiagnosticGap[!str_detect(EpiscoreDiagnosticGap$Episcore, 'Nikpay'),]


EpiscoreDiagnosticGap_Gadd$FDR <- p.adjust(EpiscoreDiagnosticGap_Gadd$P, method = 'fdr')
write.csv(EpiscoreDiagnosticGap_Gadd, file = 'dbGAP_EpiscoreGaddDiagGap30.csv')

#------------------------------------------------------------------------------#
samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]

EpiscoreDiagnosticGap <- data.frame()

for (episcore in EpiscoreNames) {
  CasesOnly[[episcore]] <- scale(CasesOnly[[episcore]])
  
  model <- glm(TimeGap ~ CasesOnly[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(TimeGap ~ CasesOnly[[episcore]], family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreDiagnosticGap <- rbind(EpiscoreDiagnosticGap, newrow)
  
}
colnames(EpiscoreDiagnosticGap) <- c('Episcore', 'Coef', 'SE', 'B', 'P')

EpiscoreDiagnosticGap_Gadd <- EpiscoreDiagnosticGap[!str_detect(EpiscoreDiagnosticGap$Episcore, 'Nikpay'),]


EpiscoreDiagnosticGap_Gadd$FDR <- p.adjust(EpiscoreDiagnosticGap_Gadd$P, method = 'fdr')
write.csv(EpiscoreDiagnosticGap_Gadd, file = 'dbGAP_EpiscoreGaddDiagGap10.csv')

#------------------------------------------------------------------------------#
samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=5,]

EpiscoreDiagnosticGap <- data.frame()

for (episcore in EpiscoreNames) {
  
  #residepiscore <- resid(glm(CasesOnly[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = CasesOnly))
  #Gap <- CasesOnly$TimeGap
  #Gap <- Gap[as.numeric(names(residepiscore))]
  CasesOnly[[episcore]] <- scale(CasesOnly[[episcore]])
  
  model <- glm(TimeGap ~ CasesOnly[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(Gap ~ residepiscore, family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreDiagnosticGap <- rbind(EpiscoreDiagnosticGap, newrow)
  
}
colnames(EpiscoreDiagnosticGap) <- c('Episcore', 'Coef', 'SE', 'B', 'P')

EpiscoreDiagnosticGap_Gadd <- EpiscoreDiagnosticGap[!str_detect(EpiscoreDiagnosticGap$Episcore, 'Nikpay'),]


EpiscoreDiagnosticGap_Gadd$FDR <- p.adjust(EpiscoreDiagnosticGap_Gadd$P, method = 'fdr')

write.csv(EpiscoreDiagnosticGap_Gadd, file = 'dbGAP_EpiscoreGaddDiagGap5.csv')

rm(EpiscoreDiagnosticGap, EpiscoreDiagnosticGap_Gadd)
#------------------------------------------------------------------------------#
#Time to diagnosis cox proportional hazards 
library(survival)
library(survminer)

samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=30,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- samples[samples$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

coxdata <- rbind(CasesOnly, ControlOnly)


EpiscoreCox_TTDiag <- data.frame()

for (episcore in EpiscoreNames) {
  

  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  res.cox <- coxph(Surv(TimeGap, CASE) ~ episcore + SEX + AGE + BMI + SMOKING, data = coxdata)

  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
   
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  
  
}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag[!str_detect(EpiscoreCox_TTDiag$Episcore, 'Nikpay'),]


EpiscoreCox_TTDiag_Gadd$FDR <- p.adjust(EpiscoreCox_TTDiag_Gadd$P, method = 'fdr')

#EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag_Gadd[-10,]

EpiscoreCox_TTDiag_Gadd[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Gadd[[x]]))
EpiscoreCox_TTDiag_Gadd[,2:7] <- signif(EpiscoreCox_TTDiag_Gadd[,2:7],3)


write.csv(EpiscoreCox_TTDiag_Gadd, file = 'dbGAP_EpiscoreGadd_Cox.csv')

PancreaticNames <- EpiscoreCox_TTDiag[EpiscoreCox_TTDiag$P<0.05,]$Episcore
rm(EpiscoreCox_TTDiag, EpiscoreCox_TTDiag_Gadd)


#------------------------------------------------------------------------------#
samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- samples[samples$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

coxdata <- rbind(CasesOnly, ControlOnly)

# Matches 
#IDMatch <- coxdata[coxdata$CASE==1]$MATCH_ID
#coxdata <- coxdata[coxdata$CASE==1 | coxdata$MATCH_ID %in% IDMatch,]

EpiscoreCox_TTDiag <- data.frame()

for (episcore in EpiscoreNames) {
  
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  res.cox <- coxph(Surv(TimeGap, CASE) ~ episcore + SEX + AGE + BMI + SMOKING, data = coxdata)
  
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  
  
}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag[!str_detect(EpiscoreCox_TTDiag$Episcore, 'Nikpay'),]


EpiscoreCox_TTDiag_Gadd$FDR <- p.adjust(EpiscoreCox_TTDiag_Gadd$P, method = 'fdr')

#EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag_Gadd[-10,]

EpiscoreCox_TTDiag_Gadd[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Gadd[[x]]))
EpiscoreCox_TTDiag_Gadd[,2:7] <- signif(EpiscoreCox_TTDiag_Gadd[,2:7],3)


write.csv(EpiscoreCox_TTDiag_Gadd, file = 'dbGAP_EpiscoreGadd_Cox10.csv')

PancreaticNames <- EpiscoreCox_TTDiag[EpiscoreCox_TTDiag$P<0.05,]$Episcore
rm(EpiscoreCox_TTDiag, EpiscoreCox_TTDiag_Gadd)

#------------------------------------------------------------------------------#
# Repeat analysis including celcounts 
samples <- samplesplaceholder
CasesOnly <- samples[samples$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- samples[samples$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

coxdata <- rbind(CasesOnly, ControlOnly)


EpiscoreCox_TTDiag <- data.frame()

for (episcore in EpiscoreNames) {
  
  
  coxdata$episcore <- coxdata[[episcore]]
  coxdata$episcore <- scale(coxdata$episcore)
  res.cox <- coxph(Surv(TimeGap, CASE) ~ episcore + SEX + AGE + BMI + SMOKING +
                     NE + EO + BA + MO + B + DC + NK, data = coxdata)
  
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  EpiscoreCox_TTDiag <- rbind(EpiscoreCox_TTDiag, newrow)
  
  
}
colnames(EpiscoreCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag[!str_detect(EpiscoreCox_TTDiag$Episcore, 'Nikpay'),]


EpiscoreCox_TTDiag_Gadd$FDR <- p.adjust(EpiscoreCox_TTDiag_Gadd$P, method = 'fdr')

#EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag_Gadd[-10,]

EpiscoreCox_TTDiag_Gadd[,2:7] <- lapply(2:7, function(x) as.numeric(EpiscoreCox_TTDiag_Gadd[[x]]))
EpiscoreCox_TTDiag_Gadd[,2:7] <- signif(EpiscoreCox_TTDiag_Gadd[,2:7],3)



#------------------------------------------------------------------------------#
#-----------------------PCClock Time to diagnosis------------------------------#
#------------------------------------------------------------------------------#
load("PCClock_dbGaP.RData")
load("dbGAPClocks.Rdata")

samples<- samples[,c('horvath', 'hannum', 'phenoage', 'brenner', 'dnamtl', 'GrimAge')]

load('dbclocks.RData')
dbclocks <- dbclocks[, c('dunedinpace', 'dunedinpoam38')]

PCClock_DNAmAge <- PCClock_DNAmAgeAccel
PCClock_DNAmAge$BMI <- ifelse(PCClock_DNAmAge$BMI<99, PCClock_DNAmAge$BMI, NA)
PCClock_DNAmAge <- cbind(PCClock_DNAmAge, samples, dbclocks)
PCClock_DNAmAge_Placeholder <- PCClock_DNAmAge
PCClock_DNAmAge_Placeholder$GrimAge <- resid(lm(PCClock_DNAmAge_Placeholder$GrimAge ~ PCClock_DNAmAge_Placeholder$AGE))
rm(PCClock_DNAmAgeAccel, samples)

PCClockDiagnosticGap <- data.frame()
ClockNames <- c('PCHorvath1Resid', 'PCHorvath2Resid', 'PCHannumResid', 'PCPhenoAgeResid', 'PCDNAmTLResid', 'PCGrimAgeResid',
                'horvath', 'hannum', 'phenoage', 'dunedinpace', 'dunedinpoam38', 'brenner', 'dnamtl', 'GrimAge')

PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder
CasesOnly <- PCClock_DNAmAge[PCClock_DNAmAge$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=30,]

for (PCClock in ClockNames) {
  
  CasesOnly[[PCClock]] <- scale(CasesOnly[[PCClock]])
  model <- glm(TimeGap ~ CasesOnly[[PCClock]] + SEX  + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(TimeGap ~ CasesOnly[[PCClock]], family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCClock, modsum)
  PCClockDiagnosticGap <- rbind(PCClockDiagnosticGap, newrow)
  CasesOnly$Clocktemp <- CasesOnly[, PCClock]
  
  #print(ggplot(CasesOnly, aes(x=TimeGap, y =Clocktemp)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  #        ggtitle(paste0('Assoc between time to diagnosis and ', PCClock)))
  
}
colnames(PCClockDiagnosticGap) <- c('PCClock', 'Coef', 'SE', 'B', 'P')
PCClockDiagnosticGap$FDR <- p.adjust(PCClockDiagnosticGap$P, method = 'fdr')
write.csv(PCClockDiagnosticGap, file = 'dbGAP_PCClock_TTD_30.csv')


#------------------------------------------------------------------------------#
PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder
CasesOnly <- PCClock_DNAmAge[PCClock_DNAmAge$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]
PCClockDiagnosticGap <- data.frame()

for (PCClock in ClockNames) {
  
  CasesOnly[[PCClock]] <- scale(CasesOnly[[PCClock]])
  
  model <- glm(TimeGap ~ CasesOnly[[PCClock]] + SEX + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(TimeGap ~ CasesOnly[[PCClock]], family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCClock, modsum)
  PCClockDiagnosticGap <- rbind(PCClockDiagnosticGap, newrow)
  CasesOnly$Clocktemp <- CasesOnly[, PCClock]
  
  #print(ggplot(CasesOnly, aes(x=TimeGap, y =Clocktemp)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  #        ggtitle(paste0('Assoc between time to diagnosis and ', PCClock)))
  
}
colnames(PCClockDiagnosticGap) <- c('PCClock', 'Coef', 'SE', 'B', 'P')
PCClockDiagnosticGap$FDR <- p.adjust(PCClockDiagnosticGap$P, method = 'fdr')
write.csv(PCClockDiagnosticGap, file = 'dbGAP_PCClock_TTD_10.csv')

#------------------------------------------------------------------------------#
PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder
CasesOnly <- PCClock_DNAmAge[PCClock_DNAmAge$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=5,]
PCClockDiagnosticGap <- data.frame()

for (PCClock in ClockNames) {
  
  CasesOnly[[PCClock]] <- scale(CasesOnly[[PCClock]])
  
  model <- glm(TimeGap ~ CasesOnly[[PCClock]] + SEX + BMI + SMOKING + RACE, family = 'gaussian', data = CasesOnly)
  #model <- glm(TimeGap ~ CasesOnly[[PCClock]], family = 'gaussian', data = CasesOnly)
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCClock, modsum)
  PCClockDiagnosticGap <- rbind(PCClockDiagnosticGap, newrow)
  CasesOnly$Clocktemp <- CasesOnly[, PCClock]
  
  #print(ggplot(CasesOnly, aes(x=TimeGap, y =Clocktemp)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  #        ggtitle(paste0('Assoc between time to diagnosis and ', PCClock)))
  
}
colnames(PCClockDiagnosticGap) <- c('PCClock', 'Coef', 'SE', 'B', 'P')
PCClockDiagnosticGap$FDR <- p.adjust(PCClockDiagnosticGap$P, method = 'fdr')
write.csv(PCClockDiagnosticGap, file = 'dbGAP_PCClock_TTD_5.csv')
rm(PCClockDiagnosticGap)
#------------------------------------------------------------------------------#
#Time to diagnosis cox proportional hazards 
library(survival)
library(survminer)

samples <- PCClock_DNAmAge_Placeholder
CasesOnly <- PCClock_DNAmAge_Placeholder[PCClock_DNAmAge_Placeholder$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=30,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- PCClock_DNAmAge_Placeholder[PCClock_DNAmAge_Placeholder$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

coxdata <- rbind(CasesOnly, ControlOnly)


PCClockCox_TTDiag <- data.frame()

for (episcore in ClockNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata[[episcore]] <- scale(coxdata[[episcore]])
  res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]] + SEX + BMI + SMOKING, data = coxdata)

  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  PCClockCox_TTDiag <- rbind(PCClockCox_TTDiag, newrow)
  
  
}
colnames(PCClockCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

PCClockCox_TTDiag$FDR <- p.adjust(PCClockCox_TTDiag$P, method = 'fdr')

#EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag_Gadd[-10,]

PCClockCox_TTDiag[,2:7] <- lapply(2:7, function(x) as.numeric(PCClockCox_TTDiag[[x]]))
PCClockCox_TTDiag[,2:7] <- signif(PCClockCox_TTDiag[,2:7],3)
write.csv(PCClockCox_TTDiag, file = 'dbGAP_PCClockTTDCox.csv')

rm(PCClockCox_TTDiag)

#------------------------------------------------------------------------------#
samples <- PCClock_DNAmAge_Placeholder
CasesOnly <- PCClock_DNAmAge_Placeholder[PCClock_DNAmAge_Placeholder$CASE==1,]
CasesOnly$TimeGap <- CasesOnly$YEAR_OF_DIAGNOSIS - CasesOnly$YEAR_OF_BLOOD_DRAWN
CasesOnly <- CasesOnly[CasesOnly$TimeGap<=10,]
recentdiag <- max(CasesOnly$YEAR_OF_DIAGNOSIS)
ControlOnly <- PCClock_DNAmAge_Placeholder[PCClock_DNAmAge_Placeholder$CASE==0,]
ControlOnly$TimeGap <- recentdiag - ControlOnly$YEAR_OF_BLOOD_DRAWN

coxdata <- rbind(CasesOnly, ControlOnly)


PCClockCox_TTDiag <- data.frame()

for (episcore in ClockNames) {
  
  #residepiscore <- resid(glm(coxdata[[episcore]] ~ AGE + SEX + BMI + SMOKING + RACE, data = coxdata))
  #coxdata[[episcore]]
  
  coxdata[[episcore]] <- scale(coxdata[[episcore]])
  res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]] + SEX + BMI + SMOKING, data = coxdata)
  
  #res.cox <- coxph(Surv(TimeGap, CASE) ~ coxdata[[episcore]], data = coxdata)
  modsum <- summary(res.cox)
  pval <- modsum[["coefficients"]][1,5]
  
  
  modsum <- modsum[["conf.int"]][1,]
  
  newrow <- c(episcore, modsum, pval)
  PCClockCox_TTDiag <- rbind(PCClockCox_TTDiag, newrow)
  
  
}
colnames(PCClockCox_TTDiag) <- c('Episcore', 'HR', 'negHR', 'LI', 'UI', 'P')

PCClockCox_TTDiag$FDR <- p.adjust(PCClockCox_TTDiag$P, method = 'fdr')

#EpiscoreCox_TTDiag_Gadd <- EpiscoreCox_TTDiag_Gadd[-10,]

PCClockCox_TTDiag[,2:7] <- lapply(2:7, function(x) as.numeric(PCClockCox_TTDiag[[x]]))
PCClockCox_TTDiag[,2:7] <- signif(PCClockCox_TTDiag[,2:7],3)
write.csv(PCClockCox_TTDiag, file = 'dbGAP_PCClockTTDCox10.csv')

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#--------------------------CASE Status associations----------------------------#
#------------------------------------------------------------------------------#
# Test basic associations between case status and episcores
EpiscoreLogReg <- data.frame()

samples <- samplesplaceholder

samples$TimeGap <- samples$YEAR_OF_DIAGNOSIS - samples$YEAR_OF_BLOOD_DRAWN 
samples <- samples[!(samples$CASE==1 & samples$TimeGap>=30),]
samples$keep <- ifelse(samples$MATCH_ID %in% samples$SUBJECT_ID.x, T, F)
samples <- samples[samples$keep==T,]
samples <- samples[samples$dbGaP_Subject_ID!='3041687']

for (episcore in EpiscoreNames) {
  
  samples[[episcore]] <- scale(samples[[episcore]])
  model <- glm(CASE ~ samples[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'binomial', data = samples)
  #model <- glm(CASE ~ samples[[episcore]], family = 'binomial', data = samples)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreLogReg <- rbind(EpiscoreLogReg, newrow)
  
}

colnames(EpiscoreLogReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')
EpiscoreLogReg_Gadd <- EpiscoreLogReg[!str_detect(EpiscoreLogReg$Episcore, 'Nikpay'),]


EpiscoreLogReg_Gadd$FDR <- p.adjust(EpiscoreLogReg_Gadd$P, method = 'fdr')
write.csv(EpiscoreLogReg_Gadd, file = 'dbGAP_EpiscoreGadd_LogReg_30.csv')

#------------------------------------------------------------------------------#
EpiscoreLogReg <- data.frame()

samples <- samplesplaceholder

samples$TimeGap <- samples$YEAR_OF_DIAGNOSIS - samples$YEAR_OF_BLOOD_DRAWN 
samples <- samples[!(samples$CASE==1 & samples$TimeGap>=10),]
samples$keep <- ifelse(samples$MATCH_ID %in% samples$SUBJECT_ID.x, T, F)
samples <- samples[samples$keep==T,]
samples <- samples[samples$dbGaP_Subject_ID!='3041687']

for (episcore in EpiscoreNames) {
  
  samples[[episcore]] <- scale(samples[[episcore]])
  
  model <- glm(CASE ~ samples[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'binomial', data = samples)
  #model <- glm(CASE ~ samples[[episcore]], family = 'binomial', data = samples)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreLogReg <- rbind(EpiscoreLogReg, newrow)
  
}

colnames(EpiscoreLogReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')
EpiscoreLogReg_Gadd <- EpiscoreLogReg[!str_detect(EpiscoreLogReg$Episcore, 'Nikpay'),]


EpiscoreLogReg_Gadd$FDR <- p.adjust(EpiscoreLogReg_Gadd$P, method = 'fdr')
write.csv(EpiscoreLogReg_Gadd, file = 'dbGAP_EpiscoreGadd_LogReg_10.csv')

#------------------------------------------------------------------------------#
EpiscoreLogReg <- data.frame()

samples <- samplesplaceholder

samples$TimeGap <- samples$YEAR_OF_DIAGNOSIS - samples$YEAR_OF_BLOOD_DRAWN 
samples <- samples[!(samples$CASE==1 & samples$TimeGap>=5),]
samples$keep <- ifelse(samples$MATCH_ID %in% samples$SUBJECT_ID.x, T, F)
samples <- samples[samples$keep==T,]
samples <- samples[samples$dbGaP_Subject_ID!='3041687']

for (episcore in EpiscoreNames) {
  
  
  samples[[episcore]] <- scale(samples[[episcore]])
  
  model <- glm(CASE ~ samples[[episcore]] + AGE + SEX + BMI + SMOKING + RACE, family = 'binomial', data = samples)
  #model <- glm(CASE ~ samples[[episcore]], family = 'binomial', data = samples)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(episcore, modsum)
  EpiscoreLogReg <- rbind(EpiscoreLogReg, newrow)
  
}

colnames(EpiscoreLogReg) <- c('Episcore', 'Coef', 'SE', 'Z', 'P')
EpiscoreLogReg_Gadd <- EpiscoreLogReg[!str_detect(EpiscoreLogReg$Episcore, 'Nikpay'),]


EpiscoreLogReg_Gadd$FDR <- p.adjust(EpiscoreLogReg_Gadd$P, method = 'fdr')
write.csv(EpiscoreLogReg_Gadd, file = 'dbGAP_EpiscoreGadd_LogReg_5.csv')
rm(EpiscoreLogReg, EpiscoreLogReg_Gadd)
#------------------------------------------------------------------------------#
#----------------------PC Clock Case association-------------------------------#
#------------------------------------------------------------------------------#

PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder

PCClock_DNAmAge$TimeGap <- PCClock_DNAmAge$YEAR_OF_DIAGNOSIS - PCClock_DNAmAge$YEAR_OF_BLOOD_DRAWN 
PCClock_DNAmAge <- PCClock_DNAmAge[!(PCClock_DNAmAge$CASE==1 & PCClock_DNAmAge$TimeGap>=30),]
PCClock_DNAmAge$keep <- ifelse(PCClock_DNAmAge$MATCH_ID %in% PCClock_DNAmAge$SUBJECT_ID.x, T, F)
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$keep==T,]
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$dbGaP_Subject_ID!='3041687',]

PCClockLogReg <- data.frame()

for (PCCLock in ClockNames) {
  
  PCClock_DNAmAge[[PCCLock]] <- scale(PCClock_DNAmAge[[PCCLock]])
  model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]] + SEX + BMI + SMOKING + RACE, family = 'binomial', data = PCClock_DNAmAge)
  #model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]], family = 'binomial', data = PCClock_DNAmAge)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCCLock, modsum)
  PCClockLogReg <- rbind(PCClockLogReg, newrow)
  
}

colnames(PCClockLogReg) <- c('PCClock', 'Coef', 'SE', 'Z', 'P')
PCClockLogReg$FDR <- p.adjust(PCClockLogReg$P, method = 'fdr')
write.csv(PCClockLogReg, file = 'dbGAP_PCClock_LogReg30.csv')

#------------------------------------------------------------------------------#
PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder

PCClock_DNAmAge$TimeGap <- PCClock_DNAmAge$YEAR_OF_DIAGNOSIS - PCClock_DNAmAge$YEAR_OF_BLOOD_DRAWN 
PCClock_DNAmAge <- PCClock_DNAmAge[!(PCClock_DNAmAge$CASE==1 & PCClock_DNAmAge$TimeGap>=10),]
PCClock_DNAmAge$keep <- ifelse(PCClock_DNAmAge$MATCH_ID %in% PCClock_DNAmAge$SUBJECT_ID.x, T, F)
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$keep==T,]
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$dbGaP_Subject_ID!='3041687',]

PCClockLogReg <- data.frame()

for (PCCLock in ClockNames) {
  
  PCClock_DNAmAge[[PCCLock]] <- scale(PCClock_DNAmAge[[PCCLock]])
  
  model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]] + SEX + AGE + BMI + SMOKING + RACE, family = 'binomial', data = PCClock_DNAmAge)
  #model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]], family = 'binomial', data = PCClock_DNAmAge)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCCLock, modsum)
  PCClockLogReg <- rbind(PCClockLogReg, newrow)
  
}

colnames(PCClockLogReg) <- c('PCClock', 'Coef', 'SE', 'Z', 'P')
PCClockLogReg$FDR <- p.adjust(PCClockLogReg$P, method = 'fdr')
write.csv(PCClockLogReg, file = 'dbGAP_PCClock_LogReg10.csv')


#------------------------------------------------------------------------------#
PCClock_DNAmAge <- PCClock_DNAmAge_Placeholder

PCClock_DNAmAge$TimeGap <- PCClock_DNAmAge$YEAR_OF_DIAGNOSIS - PCClock_DNAmAge$YEAR_OF_BLOOD_DRAWN 
PCClock_DNAmAge <- PCClock_DNAmAge[!(PCClock_DNAmAge$CASE==1 & PCClock_DNAmAge$TimeGap>=5),]
PCClock_DNAmAge$keep <- ifelse(PCClock_DNAmAge$MATCH_ID %in% PCClock_DNAmAge$SUBJECT_ID.x, T, F)
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$keep==T,]
PCClock_DNAmAge <- PCClock_DNAmAge[PCClock_DNAmAge$dbGaP_Subject_ID!='3041687',]

PCClockLogReg <- data.frame()

for (PCCLock in ClockNames) {
  
  PCClock_DNAmAge[[PCCLock]] <- scale(PCClock_DNAmAge[[PCCLock]])
  model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]] + SEX + AGE +BMI + SMOKING + RACE, family = 'binomial', data = PCClock_DNAmAge)
  #model <- glm(CASE ~ PCClock_DNAmAge[[PCCLock]], family = 'binomial', data = PCClock_DNAmAge)
  
  
  
  modsum <- summary(model)
  modsum <- modsum[["coefficients"]][2,]
  
  newrow <- c(PCCLock, modsum)
  PCClockLogReg <- rbind(PCClockLogReg, newrow)
  
}

colnames(PCClockLogReg) <- c('PCClock', 'Coef', 'SE', 'Z', 'P')
PCClockLogReg$FDR <- p.adjust(PCClockLogReg$P, method = 'fdr')
write.csv(PCClockLogReg, file = 'dbGAP_PCClock_LogReg5.csv')
rm(PCClockLogReg)

#------------------------------------------------------------------------------#
# OSCA to evaluate variance explained by episcores and clocks















