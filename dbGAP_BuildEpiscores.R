library(tidyverse)
library(data.table)
library(readxl)
library(minfi)

#------------------------------------------------------------------------------#
# Load in meth data
betas <- as.matrix(fread('Data/dbGAP.txt', sep = '\t'), rownames=1)
#betas <- as.matrix(fread('/user/work/gv20022/dbGAP.txt', sep = '\t'), rownames=1)
#betas <- as.data.frame(betas)

# Sample info
sample1 <- fread('Data/dbGAP_Sample1.txt', sep = '\t')
sample2 <- fread('Data/dbGAP_Sample2.txt', sep = '\t')
samples <- merge(sample1, sample2, by='dbGaP_Subject_ID')
rm(sample1, sample2)
samples <- as.data.frame(samples)

# order beta values to be same as sample sheet
betas <- betas[,samples$SAMPLE_ID]
betatran <- t(betas)


# Add in Age and Sex to betas
samples$Female <- ifelse(samples$SEX==2,1,0)
samples$Age <- samples$AGE
betasnames <- rownames(betas)
betasnames <- c(betasnames, 'Female', 'Age')

betas <- rbind(betas, t(samples$Female))
betas <- rbind(betas, t(samples$Age))
rownames(betas) <- betasnames

# Cell composition estimates 
library(ewastools)
CellCounts <- estimateLC(betas, ref = 'HRS')
save(CellCounts, file = 'dbGAPCellCounts.RData')

# load in PC clocks
samples$Female <- ifelse(samples$SEX==2, 1, 0)
samples$Age <- samples$AGE
# 
clocksDir <- getwd()
source(paste(clocksDir, "/run_calcPCClocks.R", sep = ""))
source(paste(clocksDir, "/run_calcPCClocks_Accel.R", sep = ""))
clocksDir <- paste(clocksDir, "/", sep = "")
PCClock_DNAmAge <- calcPCClocks(path_to_PCClocks_directory = clocksDir, datMeth = betatran, datPheno = samples)
PCClock_DNAmAgeAccel <- calcPCClocks_Accel(PCClock_DNAmAge)
save(PCClock_DNAmAge, PCClock_DNAmAgeAccel, file = 'PCClock_dbGaP.RData')

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
betanames <- rownames(betas)
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
  CGbeta <- betas[cglist,]

  # calculate episcore
  Episcore <- CGcoef*CGbeta
  Episcore <- colSums(Episcore)

  # Add to Table of data
  newrow <- c(Protein, Episcore)
  Episcores <- rbind(Episcores, newrow)
}

# Column names
Episcore_Cols <- c('Protein')
Episcore_Cols <- append(Episcore_Cols, colnames(betas))
colnames(Episcores) <- Episcore_Cols

# Make episcores gene name column the row name
Episcores <- data.frame(Episcores[,-1], row.names = Episcores[,1])
i <- c(1:ncol(betas))
Episcores[ , i] <- apply(Episcores[ , i], 2, function(x) as.numeric(as.character(x)))
Episcores <- t(Episcores)

# Bind Episcores to samplesheet
samples <- cbind(samples, Episcores)

#------------------------------------------------------------------------------#
# Add in other clocks

library(meffonym)

scores <- cbind(age=samples$Age,
                sapply(meffonym.models(), function(model) {
                  meffonym.score(as.matrix(betas), model)$score
                }))

age.scores <- scores[,c("horvath", "hannum", 'phenoage', 'brenner', 'dnamtl')]

pace.scores <- scores[,c('dunedinpace', 'dunedinpoam38')]


acc <- apply(age.scores, 2, function(estimate) {
  residuals(lm(estimate ~ scores[,"age"]))
})

dbclocks <- cbind(acc, pace.scores)
save(dbclocks, file = 'dbclocks.RData')

samples <- cbind(samples, acc, pace.scores)


#------------------------------------------------------------------------------#
# Add in Actual GRIMAGE!!
coefficients <- read.csv('GrimAgeV1.csv')

# Grab COX normaiser and intercept values
Intercepts <- coefficients[coefficients$var == 'Intercept',]
Cox <- coefficients[coefficients$Y.pred == 'COX',]


# Reduce Episcore_Coefs to CpGs in beta values
betanames <- rownames(betas)
coefficients <- coefficients[coefficients$var %in% betanames,]

# Grab relevant coefficient values
Episcore_Genes <- unique(coefficients$Y.pred)
Episcore_Genes <- Episcore_Genes[grep('DNAm', Episcore_Genes)]



# Empty table
Episcores <- data.frame()

# Loop to create Subclocks
for(Protein in Episcore_Genes){
  
  # Load in Protein specific episcore details
  df <- data.frame(coefficients[coefficients$Y.pred==Protein,])
  
  #Grab protein specific CpG sites
  cglist <- unique(df$var)
  CGcoef <- df$beta
  CGbeta <- betas[cglist,]
  
  # calculate episcore
  Episcore <- CGcoef*CGbeta
  Episcore <- colSums(Episcore)
  
  # Add to Table of data
  newrow <- c(Protein, Episcore)
  Episcores <- rbind(Episcores, newrow)
}

# Column names
Episcore_Cols <- c('Protein')
Episcore_Cols <- append(Episcore_Cols, colnames(betas))
colnames(Episcores) <- Episcore_Cols

# Make episcores gene name column the row name
Episcores <- data.frame(Episcores[,-1], row.names = Episcores[,1])
i <- c(1:ncol(betas))
Episcores[ , i] <- apply(Episcores[ , i], 2, function(x) as.numeric(as.character(x)))
Episcores <- t(Episcores)

# Add in intercept values
Episcores <- as.data.frame(Episcores)

Intercepts <- Intercepts$beta
Episcores <- Episcores %>%
  mutate(across(everything(), ~ . + Intercepts[which(names(Episcores) == cur_column())]))

# Add in Age and sex, then normalise
Episcores <- cbind(Episcores, samples$Age, samples$Female)
colnames(Episcores) <- Cox$var

CoxNorm <- Cox$beta

Episcores <- Episcores %>%
  mutate(across(everything(), ~ . * CoxNorm[which(names(Episcores) == cur_column())]))

Episcores <- Episcores %>%
  rowwise() %>%
  mutate(GrimAge = sum(c_across(everything())))


# Bind Episcores to samplesheet
samples <- cbind(samples, Episcores)


#------------------------------------------------------------------------------#

save(samples, file = 'dbGAPClocks.RData')


































