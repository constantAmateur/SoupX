## ----global_options, include=FALSE--------------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE)
knitr::opts_chunk$set(dev='png')

## ----quick_start, eval=FALSE--------------------------------------------------
#  install.packages('SoupX')
#  library(SoupX)
#  #Load data and estimate soup profile
#  sc = load10X('Path/to/cellranger/outs/folder/')
#  #Estimate rho
#  sc = autoEstCont(sc)
#  #Clean the data
#  out = adjustCounts(sc)

## ----install_CRAN,eval=FALSE--------------------------------------------------
#  install.packages('SoupX')

## ----install, eval=FALSE------------------------------------------------------
#  devtools::install_github("constantAmateur/SoupX",ref='devel')

## ----load---------------------------------------------------------------------
library(SoupX)

## ----download,eval=FALSE------------------------------------------------------
#  tmpDir = tempdir(check=TRUE)
#  download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))
#  download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))
#  untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)
#  untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)

## ----load_data,eval=FALSE-----------------------------------------------------
#  sc = load10X(tmpDir)

## ----load_data_manual,eval=FALSE----------------------------------------------
#  toc = Seurat::Read10X(file.path(tmpDir,'filtered_gene_bc_matrices','GRCh38'))
#  tod = Seurat::Read10X(file.path(tmpDir,'raw_gene_bc_matrices','GRCh38'))
#  sc = SoupChannel(tod,toc)

## ----load_saved---------------------------------------------------------------
data(PBMC_sc)
sc = PBMC_sc
sc

## ----estimateSoup, eval=FALSE-------------------------------------------------
#  sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
#  sc = estimateSoup(sc)

## ----estimateNoDrops, eval=TRUE-----------------------------------------------
library(Matrix)
toc = sc$toc
scNoDrops = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#Calculate soup profile
soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)

## ----set_clustering-----------------------------------------------------------
data(PBMC_metaData)
sc = setClusters(sc,setNames(PBMC_metaData$Cluster,rownames(PBMC_metaData)))

## ----add_DR-------------------------------------------------------------------
sc = setDR(sc,PBMC_metaData[colnames(sc$toc),c('RD1','RD2')])

## ----plot_annot---------------------------------------------------------------
library(ggplot2)
dd = PBMC_metaData[colnames(sc$toc),]
mids = aggregate(cbind(RD1,RD2) ~ Annotation,data=dd,FUN=mean)
gg = ggplot(dd,aes(RD1,RD2)) + 
  geom_point(aes(colour=Annotation),size=0.2) +
  geom_label(data=mids,aes(label=Annotation)) +
  ggtitle('PBMC 4k Annotation') +
  guides(colour = guide_legend(override.aes = list(size=1)))
plot(gg)

## ----plot_IGKC----------------------------------------------------------------
dd$IGKC = sc$toc['IGKC',]
gg = ggplot(dd,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)

## ----sanity_check-------------------------------------------------------------
gg = plotMarkerMap(sc,'IGKC')
plot(gg)

## ----set_rho------------------------------------------------------------------
sc = setContaminationFraction(sc,0.2)

## ----genes1-------------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC'))

## ----genes2-------------------------------------------------------------------
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))

## ----auto_est-----------------------------------------------------------------
sc = autoEstCont(sc)

## ----auto_est_unif_prior------------------------------------------------------
sc = autoEstCont(sc,priorRhoStdDev=0.3)

## ----topSoupGenes-------------------------------------------------------------
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)

## ----inferNonExpressed--------------------------------------------------------
plotMarkerDistribution(sc)

## ----igGenes------------------------------------------------------------------
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')

## ----calculateNullMatrix------------------------------------------------------
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes),clusters=FALSE)

## ----visNullMatrix------------------------------------------------------------
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)

## ----calcNullMatrixWithClustering---------------------------------------------
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes))
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)

## ----calcContamination--------------------------------------------------------
sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst)

## ----viewCont-----------------------------------------------------------------
head(sc$metaData)

## ----decontaminate------------------------------------------------------------
out = adjustCounts(sc)

## ----mostZeroed---------------------------------------------------------------
cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed

## ----mostReduced--------------------------------------------------------------
tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=20)

## ----IGKC_change--------------------------------------------------------------
plotChangeMap(sc,out,'IGKC')

## ----change_plots,eval=FALSE--------------------------------------------------
#  plotChangeMap(sc,out,'LYZ')
#  plotChangeMap(sc,out,'CD74')
#  plotChangeMap(sc,out,'IL32')
#  plotChangeMap(sc,out,'TRAC')
#  plotChangeMap(sc,out,'S100A9')
#  plotChangeMap(sc,out,'NKG7')
#  plotChangeMap(sc,out,'GNLY')
#  plotChangeMap(sc,out,'CD4')
#  plotChangeMap(sc,out,'CD8A')

## ----writeOut,eval=FALSE------------------------------------------------------
#  DropletUtils:::write10xCounts('./strainedCounts',out)

## ----seurat,eval=FALSE--------------------------------------------------------
#  library(Seurat)
#  srat = CreateSeuratObject(out)

## ----seuratMulti,eval=FALSE---------------------------------------------------
#  library(Seurat)
#  srat = list()
#  for(nom in names(scs)){
#    #Clean channel named 'nom'
#    tmp = adjustCounts(scs[[nom]])
#    #Add experiment name to cell barcodes to make them unique
#    colnames(tmp) = paste0(nom,'_',colnames(tmp))
#    #Store the result
#    srat[[nom]] = tmp
#  }
#  #Combine all count matricies into one matrix
#  srat = do.call(cbind,srat)
#  srat = CreateSeuratObject(srat)

