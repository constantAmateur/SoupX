## ----global_options, include=FALSE---------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE)

## ----load_scl------------------------------------------------------------
library(SoupX)
data(PBMC_DR)
data(PBMC_SCL)

## ----extractCellLabs-----------------------------------------------------
labs = rownames(PBMC_DR)[PBMC_DR$Cluster %in% c(1,2,3,5,6,8)]

## ----useToEst------------------------------------------------------------
toUse = matrix(colnames(PBMC_SCL$toc) %in% labs,nrow=1,dimnames=list('IG',colnames(PBMC_SCL$toc)))

## ----estRho--------------------------------------------------------------
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')
PBMC_SCL = calculateContaminationFraction(PBMC_SCL,'Channel1',list(IG=igGenes),useToEst=toUse)

## ----estRhoPlot1---------------------------------------------------------
plotChannelContamination(PBMC_SCL,'Channel1')

## ----estRhoPlot2---------------------------------------------------------
PBMC_SCL = calculateContaminationFraction(PBMC_SCL,'Channel1',list(IG=igGenes))
plotChannelContamination(PBMC_SCL,'Channel1')

## ----useToEstAuto--------------------------------------------------------
toUsePoisson = identifyExpressingCells(PBMC_SCL,'Channel1',list(IG=igGenes))
toUsePoisson[,1:10]
table(toUsePoisson[1,])

## ----useToEstComb--------------------------------------------------------
toUseComb = toUse & toUsePoisson
sum(toUse)
sum(toUseComb)
PBMC_SCL = calculateContaminationFraction(PBMC_SCL,'Channel1',list(IG=igGenes),useToEst = toUseComb)
plotChannelContamination(PBMC_SCL,'Channel1')

## ----smallerBins---------------------------------------------------------
PBMC_SCL = calculateContaminationFraction(PBMC_SCL,'Channel1',list(IG=igGenes),useToEst = toUseComb,tgtSoupCntsPerGroup=100)
plotChannelContamination(PBMC_SCL,'Channel1')

## ----interpolate---------------------------------------------------------
PBMC_SCL = interpolateCellContamination(PBMC_SCL,'Channel1')

## ----fixedRho------------------------------------------------------------
PBMC_SCL = interpolateCellContamination(PBMC_SCL,'Channel1',interpolationMethod='fixed',fixedContaminationValue=0.05)

## ----inferGenes----------------------------------------------------------
PBMC_SCL = inferNonExpressedGenes(PBMC_SCL)
plotMarkerDistribution(PBMC_SCL,'Channel1')

