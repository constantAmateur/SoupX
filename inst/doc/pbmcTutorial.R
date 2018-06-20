## ----global_options, include=FALSE---------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE)

## ----genes1--------------------------------------------------------------
nonExpressedGeneList = list(HEM=c('HBB','HBA2'),IG = c('IGKC'))

## ----genes2--------------------------------------------------------------
nonExpressedGeneList = list(HEM=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))

## ----install, eval=FALSE-------------------------------------------------
#  devtools::install_github("constantAmateur/SoupX")

## ----load----------------------------------------------------------------
library(SoupX)

## ----load_data-----------------------------------------------------------
library(SoupX)
dataDirs = c('SoupX_pbmc4k_demo/')
scl = load10X(dataDirs)

## ----estimateSoup, eval=FALSE--------------------------------------------
#  scl = load10X(dataDirs,keepDroplets=TRUE)
#  scl$channels$Channel1 = estimateSoup(scl$channels$Channel1)

## ----init_dataset--------------------------------------------------------
data(PBMC_DR)

## ----plot_IGKC-----------------------------------------------------------
library(ggplot2)
PBMC_DR$IGKC = scl$toc['IGKC',rownames(PBMC_DR)]
gg = ggplot(PBMC_DR,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)

## ----sanity_check--------------------------------------------------------
gg = plotMarkerMap(scl,'IGKC',PBMC_DR)
plot(gg)

## ----inferNonExpressed---------------------------------------------------
scl = inferNonExpressedGenes(scl)

## ----plotNonExpressed,fig.width=12---------------------------------------
tstGenes = rownames(scl$channels$Channel1$nonExpressedGenes)[seq(20)]
gg = plotMarkerDistribution(scl,'Channel1',tstGenes)
plot(gg)

## ----plotAlt-------------------------------------------------------------
gg = plotMarkerDistribution(scl,'Channel1')

## ----igGenes-------------------------------------------------------------
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')

## ----calcContamination---------------------------------------------------
scl = calculateContaminationFraction(scl,'Channel1',list(IG=igGenes))
gg = plotChannelContamination(scl,'Channel1')
plot(gg)

## ----hgGenes-------------------------------------------------------------
hgGenes = c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
scl = calculateContaminationFraction(scl,'Channel1',list(IG=igGenes,HG=hgGenes))

## ----hgNULL--------------------------------------------------------------
library(Matrix)
rowSums(scl$channels$Channel1$toc[hgGenes,])

## ----interpolate---------------------------------------------------------
scl = interpolateCellContamination(scl,'Channel1',useGlobal=TRUE)

## ----rhoHead-------------------------------------------------------------
head(scl$channels$Channel1$rhos)

## ----doCorrect-----------------------------------------------------------
scl = strainCells(scl)
scl = adjustCounts(scl)

## ----mostZeroed----------------------------------------------------------
cntSoggy = rowSums(scl$toc>0)
cntStrained = rowSums(scl$strainedExp>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed

## ----mostAdjusted--------------------------------------------------------
cntAdjusted = rowSums(scl$atoc>0)
((cntSoggy-cntAdjusted)/cntSoggy)[names(mostZeroed)]

## ----mostStrained--------------------------------------------------------
#Need to convert table of counts to a normalised expression matrix
soggyExp = t(t(scl$toc)/scl$nUMIs)
rowSums(scl$strainedExp[names(mostZeroed),] < soggyExp[names(mostZeroed),])/rowSums(scl$toc[names(mostZeroed),]>0)

## ----mostDecreased-------------------------------------------------------
tail(sort(rowSums(scl$strainedExp < soggyExp)),n=10)

## ----MTchange------------------------------------------------------------
(cntSoggy-cntAdjusted)[grep('^MT-',names(cntAdjusted))]

## ----IGKC_change, fig.height=8, fig.width=20-----------------------------
plotChangeMap(scl,'IGKC',PBMC_DR)

## ----IGKC_change_ratio, fig.height=8, fig.width=20-----------------------
plotChangeMap(scl,'IGKC',PBMC_DR,dataType='ratio')

## ----change_plots, fig.height=8, fig.width=16----------------------------
plotChangeMap(scl,'LYZ',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'CD74',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'HLA-DRA',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'IL32',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'TRAC',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'CD3D',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'S100A9',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'S100A8',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'LTB',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'NKG7',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'GNLY',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'CD4',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
plotChangeMap(scl,'CD8A',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))

## ----Seurat--------------------------------------------------------------
srat = createCleanedSeurat(scl)
srat

## ----documentation_figures, include=FALSE--------------------------------
#Some extra bits of code for making example plots
#Basic annotation
cMap = c('0'='MNP',
         '1'='CD8 T-Cell',
         '2'='CD8 T-Cell',
         '3'='CD4 T-Cell',
         '4'='B-Cell',
         '5'='CD4 T-Cell',
         '6'='NK',
         '7'='B-Cell',
         '8'='NK',
         '9'='MNP',
         '10'='MNP',
         '11'='?')
PBMC_DR$Annotation = factor(cMap[as.character(PBMC_DR$Cluster)])
mids = lapply(split(PBMC_DR[,1:2],PBMC_DR$Annotation),apply,2,mean)
mids = cbind(as.data.frame(do.call(rbind,mids)),Annotation=names(mids))
mids[1,1:2]= mids[1,1:2]+5
gg = ggplot(PBMC_DR,aes(RD1,RD2)) +
  geom_point(aes(colour=Annotation)) +
  geom_label(data=mids,aes(label=Annotation),size=16)+
  guides(colour=FALSE)+
  theme_grey(base_size = 36) +
  xlab('tSNE1')+
  ylab('tSNE2')
#This assumes this is being run with rmarkdown::render with defaults, which places the current working directory to the vignette directory.
png('../inst/images/PBMC_Annotation.png',width=960,height=960)
plot(gg)
dev.off()
#Tarted up before and after shots
gg = plotChangeMap(scl,'IGKC',PBMC_DR,includePanels=c('Uncorrected','CorrectedCounts'))
gg = ggplot(gg$df,aes(RD1,RD2)) +
  geom_point(data=gg$df[!gg$df$data,],colour='#808080',size=2) +
  geom_point(data=gg$df[gg$df$data,],colour='#e60000',size=4) +
  facet_grid(~correction) +
  xlab('tSNE1') +
  ylab('tSNE2') +
  theme_grey(base_size = 36) +
  ggtitle("IGKC expression (red) before and after decontamination")
png('../inst/images/IGKC_comparison.png',width=1920,height=960)
plot(gg)
dev.off()

