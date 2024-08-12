#EPIC quick meffonym script
library(tidyverse)
library(readr)

options(mc.cores=6)

#------------------------------------------------------------------------------#

load('EPIC_Italy_Meffil.RData')
load('EPIC_Italy_Samples.RData')

beta <- beta$beta
samplesheet <- samplesheet[samplesheet$Sample_Name %in% colnames(beta),]

# Build PC clocks
samplesheet$Female <- ifelse(samplesheet$Sex==' F', 1, 0)
samplesheet$Age <- as.numeric(samplesheet$Age)


library(meffonym)

scores <- cbind(age=samplesheet$Age,
                sapply(meffonym.models(), function(model) {
                  meffonym.score(beta, model)$score
                }))

age.scores <- scores[,c("horvath", "hannum", 'phenoage', 'brenner', 'dnamtl')]

pace.scores <- scores[,c('dunedinpace', 'dunedinpoam38')]


acc <- apply(age.scores, 2, function(estimate) {
  residuals(lm(estimate ~ scores[,"age"]))
})

acc <- cbind(acc, pace.scores)

save(acc, file = 'EPICItaly_StandardClocks.RData')


# Cell composition estimates 
library(ewastools)
CellCounts <- estimateLC(beta, ref = 'HRS')
save(CellCounts, file = 'EPICItalyCellCounts.RData')