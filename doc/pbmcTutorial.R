## ----global_options, include=FALSE---------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE)

## ----quick_start, eval=FALSE---------------------------------------------
#  library(SoupX)
#  #Load data and estimate soup profile
#  sc = load10X('Path/to/cellranger/outs/folder/')
#  #Estimate rho
#  sc = autoEstCont(sc)
#  #Clean teh data
#  out = adjustCounts(sc)

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

## ----set_clustering------------------------------------------------------
sc = setClusters(sc,PBMC_DR$Cluster)

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

## ----add_DR--------------------------------------------------------------
sc = setDR(sc,PBMC_DR)

## ----sanity_check_auto---------------------------------------------------
gg = plotMarkerMap(sc,'IGKC')
plot(gg)

## ----set_rho-------------------------------------------------------------
sc = setContaminationFraction(sc,0.2)

## ----genes1--------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC'))

## ----genes2--------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))

## ----auto_est------------------------------------------------------------
sc = autoEstCont(sc)

## ----auto_est_unif_prior-------------------------------------------------
sc = autoEstCont(sc,priorRhoStdDev=0.3)

## ----auto_est_stupid_prior-----------------------------------------------
sc = autoEstCont(sc,priorRho=0.3,priorRhoStdDev=0.05)

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

## ----decontaminate-------------------------------------------------------
out = adjustCounts(sc,clusters = setNames(PBMC_DR$Cluster,rownames(PBMC_DR)))

## ----mostZeroed----------------------------------------------------------
library(Matrix)
cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed

## ----mostReduced---------------------------------------------------------
tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=20)

## ----IGKC_change---------------------------------------------------------
plotChangeMap(sc,out,'IGKC')

## ----change_plots--------------------------------------------------------
plotChangeMap(sc,out,'LYZ')
plotChangeMap(sc,out,'CD74')
plotChangeMap(sc,out,'HLA-DRA')
plotChangeMap(sc,out,'IL32')
plotChangeMap(sc,out,'TRAC')
plotChangeMap(sc,out,'CD3D')
plotChangeMap(sc,out,'S100A9')
plotChangeMap(sc,out,'S100A8')
plotChangeMap(sc,out,'LTB')
plotChangeMap(sc,out,'NKG7')
plotChangeMap(sc,out,'GNLY')
plotChangeMap(sc,out,'CD4')
plotChangeMap(sc,out,'CD8A')

## ----writeOut------------------------------------------------------------
DropletUtils:::write10xCounts('./strainedCounts',out)

