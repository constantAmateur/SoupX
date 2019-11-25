## ----global_options, include=FALSE---------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE)

## ----genes1--------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC'))

## ----genes2--------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))

## ----install, eval=FALSE-------------------------------------------------
#  devtools::install_github("constantAmateur/SoupX")

## ----load----------------------------------------------------------------
library(SoupX)

## ----load_data-----------------------------------------------------------
library(SoupX)
dataDirs = c('SoupX_pbmc4k_demo/')
sc = load10X(dataDirs)

## ----estimateSoup, eval=FALSE--------------------------------------------
#  sc = load10X(dataDirs,keepDroplets=TRUE)
#  sc = estimateSoup(sc)

## ----init_dataset--------------------------------------------------------
data(PBMC_DR)

## ----plot_IGKC-----------------------------------------------------------
library(ggplot2)
PBMC_DR$IGKC = sc$toc['IGKC',rownames(PBMC_DR)]
gg = ggplot(PBMC_DR,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)

## ----sanity_check--------------------------------------------------------
gg = plotMarkerMap(sc,'IGKC',PBMC_DR)
plot(gg)

## ----topSoupGenes--------------------------------------------------------
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)

## ----inferNonExpressed---------------------------------------------------
plotMarkerDistribution(sc)

## ----igGenes-------------------------------------------------------------
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')

## ----calculateNullMatrix-------------------------------------------------
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes))

## ----visNullMatrix-------------------------------------------------------
plotMarkerMap(sc,geneSet=igGenes,DR=PBMC_DR,useToEst=useToEst)

## ----calcNullMatrixWithClustering----------------------------------------
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes),clusters=setNames(PBMC_DR$Cluster,rownames(PBMC_DR)))
plotMarkerMap(sc,geneSet=igGenes,DR=PBMC_DR,useToEst=useToEst)

## ----calcContamination---------------------------------------------------
sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst)

## ----viewCont------------------------------------------------------------
head(sc$metaData)

## ----cellSpecificRho,eval=FALSE------------------------------------------
#  sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst,cellSpecificEstimates=TRUE)
#  quantile(sc$metaData$rho)

## ----manualRho, eval=FALSE-----------------------------------------------
#  sc = setContaminationFraction(sc,0.1)

## ----decontaminate-------------------------------------------------------
out = adjustCounts(sc)

## ----mostZeroed----------------------------------------------------------
library(Matrix)
cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed

## ----mostReduced---------------------------------------------------------
tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=20)

## ----IGKC_change---------------------------------------------------------
plotChangeMap(sc,out,'IGKC',PBMC_DR)

## ----change_plots--------------------------------------------------------
plotChangeMap(sc,out,'LYZ',PBMC_DR)
plotChangeMap(sc,out,'CD74',PBMC_DR)
plotChangeMap(sc,out,'HLA-DRA',PBMC_DR)
plotChangeMap(sc,out,'IL32',PBMC_DR)
plotChangeMap(sc,out,'TRAC',PBMC_DR)
plotChangeMap(sc,out,'CD3D',PBMC_DR)
plotChangeMap(sc,out,'S100A9',PBMC_DR)
plotChangeMap(sc,out,'S100A8',PBMC_DR)
plotChangeMap(sc,out,'LTB',PBMC_DR)
plotChangeMap(sc,out,'NKG7',PBMC_DR)
plotChangeMap(sc,out,'GNLY',PBMC_DR)
plotChangeMap(sc,out,'CD4',PBMC_DR)
plotChangeMap(sc,out,'CD8A',PBMC_DR)

## ----writeOut------------------------------------------------------------
DropletUtils:::write10xCounts('./strainedCounts',out)

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
gg = plotChangeMap(sc,out,'IGKC',PBMC_DR)
gg = gg +
  xlab('tSNE1') +
  ylab('tSNE2') +
  ggtitle("log10 IGKC expression change after decontamination")
png('../inst/images/IGKC_comparison.png',width=1920,height=960)
plot(gg)
dev.off()

