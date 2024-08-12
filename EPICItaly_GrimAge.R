
library(tidyverse)
library(dplyr)

load('EPIC_Italy_Meffil.RData')
load('EPIC_Italy_Samples.RData')

betas <- beta$beta
samplesheet <- samplesheet[samplesheet$Sample_Name %in% colnames(betas),]
samplesheet$sex <- ifelse(samplesheet$Sex == ' F', 1, 0)
samplesheet$AGE <- samplesheet$Age
samplesheet$Age <- as.numeric(samplesheet$Age)
samplesheet$Female <- samplesheet$sex
samplesheet$SampleID <- samplesheet$Sample_Name
rm(beta)

betasnames <- rownames(betas)
betasnames <- c(betasnames, 'Female', 'Age')

betas <- rbind(betas, t(samplesheet$Female))
betas <- rbind(betas, t(samplesheet$Age))
rownames(betas) <- betasnames


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
Episcores <- cbind(Episcores, samplesheet$Age, samplesheet$Female)
colnames(Episcores) <- Cox$var

CoxNorm <- Cox$beta

Episcores <- Episcores %>%
  mutate(across(everything(), ~ . * CoxNorm[which(names(Episcores) == cur_column())]))

# Sum up composite score
Episcores <- Episcores %>%
  rowwise() %>%
  mutate(GrimAge = sum(c_across(everything())))


samplesheet <- cbind(samplesheet, Episcores)

save(samplesheet, file = 'EPIC_Italy_GrimAge.RData')




