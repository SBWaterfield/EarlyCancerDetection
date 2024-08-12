library(ggplot2)
library(readr)
library(tidyverse)


#==============================================================================#

# Clocks Cox
dbGAP_PCClockCox <- read_csv("dbGAP_PCClockTTDCox.csv")
dbGAP_PCClockCox10 <- read_csv("dbGAP_PCClockTTDCox10.csv")
EPICBreast_PCClock_Cox <- read_csv("EPICBreast_EpiscoreClock_Cox.csv")
EPICCCR_PCClock_Cox <- read_csv("EPICCCR_EpiscoreClock_Cox.csv")


Clocks <- c('horvath', 'hannum','phenoage', 'dunedinpace', 'dunedinpoam38', 'brenner', 'dnamtl', 'GrimAge',
            'horvathResid', 'hannumResid', 'phenoageResid', 'dunedinpaceResid', 'dunedinpoam38Resid', 'brennerResid', 'dnamtlResid', 'GrimAgeResid')

dbGAP_PCClockCox <- dbGAP_PCClockCox[dbGAP_PCClockCox$Episcore %in% Clocks,]
dbGAP_PCClockCox$FDR <- p.adjust(dbGAP_PCClockCox$FDR, 'fdr')

dbGAP_PCClockCox10 <- dbGAP_PCClockCox10[dbGAP_PCClockCox10$Episcore %in% Clocks,]
dbGAP_PCClockCox10$FDR <- p.adjust(dbGAP_PCClockCox10$FDR, 'fdr')

EPICBreast_PCClock_Cox <- EPICBreast_PCClock_Cox[EPICBreast_PCClock_Cox$Episcore %in% Clocks,]
EPICBreast_PCClock_Cox$FDR <- p.adjust(EPICBreast_PCClock_Cox$FDR, 'fdr')

EPICCCR_PCClock_Cox <- EPICCCR_PCClock_Cox[EPICCCR_PCClock_Cox$Episcore %in% Clocks,]
EPICCCR_PCClock_Cox$FDR <- p.adjust(EPICCCR_PCClock_Cox$FDR, 'fdr')


dbGAP_PCClockCox$Episcore <- gsub('Resid', '', dbGAP_PCClockCox$Episcore)
dbGAP_PCClockCox10$Episcore <- gsub('Resid', '', dbGAP_PCClockCox10$Episcore)
EPICBreast_PCClock_Cox$Episcore <- gsub('Resid', '', EPICBreast_PCClock_Cox$Episcore)
EPICCCR_PCClock_Cox$Episcore <- gsub('Resid', '', EPICCCR_PCClock_Cox$Episcore)

dbGAP_PCClockCox$Group <- 'Pancreatic'
dbGAP_PCClockCox10$Group <- 'Pancreatic (TTD <10 Years)'
EPICBreast_PCClock_Cox$Group <- 'Breast'
EPICCCR_PCClock_Cox$Group <- 'Colorectal'

Finaldf <- rbind(dbGAP_PCClockCox, dbGAP_PCClockCox10, EPICBreast_PCClock_Cox, EPICCCR_PCClock_Cox)
# Plot data



ggplot(Finaldf, aes(x=HR, y=Group, col=Episcore,fill=Episcore)) + geom_vline(xintercept=1) +
         geom_errorbar(aes(xmin = LI, xmax = UI), width = 0.25, position = position_dodge(width = 0.4)) +
         geom_point(position = position_dodge(width = 0.4)) + 
         ggtitle('Age Acceleration and Cancer Risk')  + xlab('HR (95% CI)') + ylab('Cancer') + 
         theme(plot.title = element_text(hjust = 0.5))

#==============================================================================#

# Episcores Cox

# Episcore results
dbGAP_EpiscoreGadd_Cox <- read_csv("dbGAP_EpiscoreGadd_Cox.csv")
dbGAP_EpiscoreGadd_Cox <- dbGAP_EpiscoreGadd_Cox[dbGAP_EpiscoreGadd_Cox$Episcore %in% c('IGFBP1', 'STC1', 'SHBG', 'SELE', 'ADIPOQ', 'ACY1'),]
#dbGAP_EpiscoreGadd_Cox <- dbGAP_EpiscoreGadd_Cox[dbGAP_EpiscoreGadd_Cox$Episcore %in% c('IGFBP1'),]

dbGAP_EpiscoreGadd_Cox$Group <- 'Pancreatic'

dbGAP_EpiscoreGadd_Cox10 <- read_csv("dbGAP_EpiscoreGadd_Cox10.csv")
dbGAP_EpiscoreGadd_Cox10 <- dbGAP_EpiscoreGadd_Cox10[dbGAP_EpiscoreGadd_Cox10$FDR <0.05,]
dbGAP_EpiscoreGadd_Cox10$Group <- 'Pancreatic (TTD <10 Years)'


EPICCR_EpiscoreGadd_Cox <- read_csv("EPICCCR_EpiscoreGadd_Cox.csv")
EPICCR_EpiscoreGadd_Cox <- EPICCR_EpiscoreGadd_Cox[EPICCR_EpiscoreGadd_Cox$Episcore %in% c('NCAM1'),]
EPICCR_EpiscoreGadd_Cox$Group <- 'Colorectal'

EPIBreast_EpiscoreGadd_Cox <- read_csv("EPICBreast_EpiscoreGadd_Cox.csv")
EPIBreast_EpiscoreGadd_Cox <- EPIBreast_EpiscoreGadd_Cox[EPIBreast_EpiscoreGadd_Cox$Episcore %in% c('INSR'),]
EPIBreast_EpiscoreGadd_Cox$Group <- 'Breast'


plotdf <- rbind(EPIBreast_EpiscoreGadd_Cox, EPICCR_EpiscoreGadd_Cox, dbGAP_EpiscoreGadd_Cox, dbGAP_EpiscoreGadd_Cox10) 
#plotdf <- rbind(EPIBreast_EpiscoreGadd_Cox, EPICCR_EpiscoreGadd_Cox, dbGAP_EpiscoreGadd_Cox) 


library(ggplot2)
ggplot(plotdf, aes(x=HR, y=Episcore, col=Group,fill=Group)) + geom_vline(xintercept=1) +
  geom_errorbar(aes(xmin = LI, xmax = UI), width = 0.25, position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) + 
  ggtitle('DNAm Protein Estimates and Cancer Risk')  + xlab('HR (95% CI)') + ylab('Protein') + 
  theme(plot.title = element_text(hjust = 0.5))

SELE2SMR <- c('SELE', exp(-0.1285), exp(-0.1285) - 1.96 * 0.036, exp(-0.1285) + 1.96 * 0.036, 'Pancreatic_2SMR')
SELEPlot <- plotdf[plotdf$Episcore=='SELE',]
SELEPlot <- SELEPlot[c('Episcore', 'HR', 'LI', 'UI', 'Group')]
SELEPlot <- rbind(SELEPlot, SELE2SMR)
SELEPlot[c(2:4)] <- lapply(SELEPlot[c(2:4)], as.numeric)

ggplot(SELEPlot, aes(x=HR, y=Episcore, col=Group,fill=Group)) + geom_vline(xintercept=1) +
  geom_errorbar(aes(xmin = LI, xmax = UI), width = 0.25, position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) + 
  ggtitle('DNAm Protein Estimates and Cancer Risk')  + xlab('HR (95% CI)') + ylab('Protein') + 
  theme(plot.title = element_text(hjust = 0.5))
#==============================================================================#
#==============================================================================#
# Time to diagnosis in clocks

library(readr)
dbGAP_PCClock_TTD_30 <- read_csv("dbGAP_PCClock_TTD_30.csv")
dbGAP_PCClock_TTD_10 <- read_csv("dbGAP_PCClock_TTD_10.csv")
EPICBreast_PCClock_TTD <- read_csv("EPICBreast_EpiscoreClock_TTD.csv")
EPICCCR_PCClock_TTD <- read_csv("EPICCCR_EpiscoreClock_TTD.csv")

dbGAP_PCClock_TTD_30 <- dbGAP_PCClock_TTD_30[dbGAP_PCClock_TTD_30$PCClock %in% Clocks,]
dbGAP_PCClock_TTD_10 <- dbGAP_PCClock_TTD_10[dbGAP_PCClock_TTD_10$PCClock %in% Clocks,]
EPICBreast_PCClock_TTD <- EPICBreast_PCClock_TTD[EPICBreast_PCClock_TTD$Episcore %in% Clocks,]
EPICCCR_PCClock_TTD <- EPICCCR_PCClock_TTD[EPICCCR_PCClock_TTD$Episcore %in% Clocks,]

dbGAP_PCClock_TTD_30$FDR <- p.adjust(dbGAP_PCClock_TTD_30$P, 'fdr')
dbGAP_PCClock_TTD_10$FDR <- p.adjust(dbGAP_PCClock_TTD_10$P, 'fdr')
EPICBreast_PCClock_TTD$FDR <- p.adjust(EPICBreast_PCClock_TTD$P, 'fdr')
EPICCCR_PCClock_TTD$FDR <- p.adjust(EPICCCR_PCClock_TTD$P, 'fdr')


dbGAP_PCClock_TTD_30$LI <- dbGAP_PCClock_TTD_30$Coef - dbGAP_PCClock_TTD_30$SE * 1.96
dbGAP_PCClock_TTD_30$UI <- dbGAP_PCClock_TTD_30$Coef + dbGAP_PCClock_TTD_30$SE * 1.96
dbGAP_PCClock_TTD_30$Group <- 'Pancreatic'
colnames(dbGAP_PCClock_TTD_30) <- c('row', 'Clock', 'Coef', 'SE', 'B', 'P', 'FDR', 'LI', 'UI', 'Group')
dbGAP_PCClock_TTD_30 <- dbGAP_PCClock_TTD_30[, c('Clock', 'Coef', 'SE',  'P', 'LI', 'UI', 'FDR', 'Group')]

dbGAP_PCClock_TTD_10$LI <- dbGAP_PCClock_TTD_10$Coef - dbGAP_PCClock_TTD_10$SE * 1.96
dbGAP_PCClock_TTD_10$UI <- dbGAP_PCClock_TTD_10$Coef + dbGAP_PCClock_TTD_10$SE * 1.96
dbGAP_PCClock_TTD_10$Group <- 'Pancreatic (TTD <10 Years)'
colnames(dbGAP_PCClock_TTD_10) <- c('row', 'Clock', 'Coef', 'SE', 'B', 'P', 'FDR', 'LI', 'UI', 'Group')
dbGAP_PCClock_TTD_10 <- dbGAP_PCClock_TTD_10[, c('Clock', 'Coef', 'SE',  'P', 'LI', 'UI', 'FDR', 'Group')]

EPICBreast_PCClock_TTD$LI <- EPICBreast_PCClock_TTD$Coef - EPICBreast_PCClock_TTD$SE * 1.96
EPICBreast_PCClock_TTD$UI <- EPICBreast_PCClock_TTD$Coef + EPICBreast_PCClock_TTD$SE * 1.96
EPICBreast_PCClock_TTD$Group <- 'Breast'
colnames(EPICBreast_PCClock_TTD) <- c('row', 'Clock', 'Coef', 'SE', 'B', 'P', 'OR', 'LI', 'UI', 'FDR', 'Group')
EPICBreast_PCClock_TTD <- EPICBreast_PCClock_TTD[, c('Clock', 'Coef', 'SE',  'P', 'LI', 'UI', 'FDR', 'Group')]
EPICBreast_PCClock_TTD$Clock <- gsub('Resid', '', EPICBreast_PCClock_TTD$Clock)

EPICCCR_PCClock_TTD$LI <- EPICCCR_PCClock_TTD$Coef - EPICCCR_PCClock_TTD$SE * 1.96
EPICCCR_PCClock_TTD$UI <- EPICCCR_PCClock_TTD$Coef + EPICCCR_PCClock_TTD$SE * 1.96
EPICCCR_PCClock_TTD$Group <- 'CRC'
colnames(EPICCCR_PCClock_TTD) <- c('row', 'Clock', 'Coef', 'SE', 'B', 'P', 'FDR', 'LI', 'UI', 'Group')
EPICCCR_PCClock_TTD <- EPICCCR_PCClock_TTD[, c('Clock', 'Coef', 'SE',  'P', 'LI', 'UI', 'FDR', 'Group')]
EPICCCR_PCClock_TTD$Clock <- gsub('Resid', '', EPICCCR_PCClock_TTD$Clock)


Finaldf <- rbind(dbGAP_PCClock_TTD_30, dbGAP_PCClock_TTD_10, EPICBreast_PCClock_TTD, EPICCCR_PCClock_TTD)
ggplot(Finaldf, aes(x=Coef, y=Group, col=Clock,fill=Clock)) + geom_vline(xintercept=0) +
  geom_errorbar(aes(xmin = LI, xmax = UI), width = 0.25, position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) + 
  ggtitle('Age Acceleration and Time to diagnosis')  + xlab('Coef (95% CI)') + ylab('Cancer') + 
  theme(plot.title = element_text(hjust = 0.5))
#==============================================================================#

# TTD in episcores

dbGAP_EpiscoreGaddDiagGap30 <- read_csv("dbGAP_EpiscoreGaddDiagGap30.csv")
dbGAP_EpiscoreGaddDiagGap10 <- read_csv("dbGAP_EpiscoreGaddDiagGap10.csv")
EPICCR_EpiscoreGadd_TTD <- read_csv("EPICCCR_EpiscoreGadd_TTD.csv")
EPICBreast_EpiscoreGadd_TTD <- read_csv("EPICBreast_EpiscoreGadd_TTD.csv")


EPICBreast_EpiscoreGadd_TTD <- EPICBreast_EpiscoreGadd_TTD[EPICBreast_EpiscoreGadd_TTD$FDR <0.05,]
EPICBreast_EpiscoreGadd_TTD$LI <- EPICBreast_EpiscoreGadd_TTD$Coef - EPICBreast_EpiscoreGadd_TTD$SE * 1.96
EPICBreast_EpiscoreGadd_TTD$UI <- EPICBreast_EpiscoreGadd_TTD$Coef + EPICBreast_EpiscoreGadd_TTD$SE * 1.96
EPICBreast_EpiscoreGadd_TTD$Group <- 'Breast'
EPICBreast_EpiscoreGadd_TTD <- EPICBreast_EpiscoreGadd_TTD %>% select(- OR)

dbGaP30 <- dbGAP_EpiscoreGaddDiagGap30[dbGAP_EpiscoreGaddDiagGap30$FDR <0.05,]$Episcore
dbGaP10 <- dbGAP_EpiscoreGaddDiagGap10[dbGAP_EpiscoreGaddDiagGap10$FDR <0.05,]$Episcore
dbGaPNames <- unique(c(dbGaP30, dbGaP10))


dbGAP_EpiscoreGaddDiagGap10 <- dbGAP_EpiscoreGaddDiagGap10[dbGAP_EpiscoreGaddDiagGap10$Episcore %in% dbGaPNames,]
dbGAP_EpiscoreGaddDiagGap10$LI <- dbGAP_EpiscoreGaddDiagGap10$Coef - dbGAP_EpiscoreGaddDiagGap10$SE * 1.96
dbGAP_EpiscoreGaddDiagGap10$UI <- dbGAP_EpiscoreGaddDiagGap10$Coef + dbGAP_EpiscoreGaddDiagGap10$SE * 1.96
dbGAP_EpiscoreGaddDiagGap10$Group <- 'Pancreatic (TTD <10 Years)'


dbGAP_EpiscoreGaddDiagGap30 <- dbGAP_EpiscoreGaddDiagGap30[dbGAP_EpiscoreGaddDiagGap30$Episcore %in% dbGaPNames,]
dbGAP_EpiscoreGaddDiagGap30$LI <- dbGAP_EpiscoreGaddDiagGap30$Coef - dbGAP_EpiscoreGaddDiagGap30$SE * 1.96
dbGAP_EpiscoreGaddDiagGap30$UI <- dbGAP_EpiscoreGaddDiagGap30$Coef + dbGAP_EpiscoreGaddDiagGap30$SE * 1.96
dbGAP_EpiscoreGaddDiagGap30$Group <- 'Pancreatic'

colnames(dbGAP_EpiscoreGaddDiagGap10)
colnames(EPICBreast_EpiscoreGadd_TTD) <- c("...1","Episcore","Coef","SE","B","P","LI","UI", 'FDR', "Group")

plotdf <- rbind(EPICBreast_EpiscoreGadd_TTD, dbGAP_EpiscoreGaddDiagGap30, dbGAP_EpiscoreGaddDiagGap10)


ggplot(plotdf, aes(x=Coef, y=Episcore, col=Group,fill=Group)) + geom_vline(xintercept=0) +
  geom_errorbar(aes(xmin = LI, xmax = UI), width = 0.25, position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4)) + 
  ggtitle('DNAm Protein Estimates and Time to Diagnosis')  + xlab('Coef (95% CI)') + ylab('Protein') + 
  theme(plot.title = element_text(hjust = 0.5))


































