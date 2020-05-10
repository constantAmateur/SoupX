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

## ----download,eval=FALSE-------------------------------------------------
#  tmpDir = tempdir(check=TRUE)
#  download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))
#  download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))
#  untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)
#  untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)

## ----load_data,eval=FALSE------------------------------------------------
#  sc = load10X(tmpDir)

## ----load_data_manual,eval=TRUE------------------------------------------
library(Matrix)
#Get the sparse matrix
tod = gzcon(url('https://github.com/constantAmateur/SoupX/raw/devel/inst/extdata/PBMC_4k/raw_gene_bc_matrices/GRCh38/matrix.mtx.gz'))
tod = readMM(tod)
#Get the gene names
tmp = gzcon(url('https://github.com/constantAmateur/SoupX/raw/devel/inst/extdata/PBMC_4k/raw_gene_bc_matrices/GRCh38/genes.tsv.gz'))
tmp = readLines(tmp)
rownames(tod) = make.unique(do.call(rbind,strsplit(tmp,'\t'))[,2])
#And the droplet barcodes
tmp = gzcon(url('https://github.com/constantAmateur/SoupX/raw/devel/inst/extdata/PBMC_4k/raw_gene_bc_matrices/GRCh38/barcodes.tsv.gz'))
colnames(tod) = gsub('-1$','',readLines(tmp))
#Load PBMC data to get cellular barcodes
toc = tod[,PBMC_cellBarcodes]
sc = SoupChannel(tod,toc)

## ----estimateSoup, eval=TRUE---------------------------------------------
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
sc = estimateSoup(sc)

## ----estimateNoDrops, eval=TRUE------------------------------------------
scNoDrops = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#Calculate soup profile
soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)

## ----load_data_with_mDat, eval=TRUE--------------------------------------
scWithMdat = SoupChannel(tod,toc,metaData=PBMC_DR)

## ----set_clustering------------------------------------------------------
sc = setClusters(sc,PBMC_DR$Cluster)

## ----add_DR--------------------------------------------------------------
sc = setDR(sc,PBMC_DR)

## ----plot_annot----------------------------------------------------------
library(ggplot2)
mids = aggregate(cbind(RD1,RD2) ~ Annotation,data=PBMC_DR,FUN=mean)
gg = ggplot(PBMC_DR,aes(RD1,RD2)) + 
  geom_point(aes(colour=Annotation),size=0.2) +
  geom_label(data=mids,aes(label=Annotation)) +
  ggtitle('PBMC 4k Annotation') +
  guides(colour = guide_legend(override.aes = list(size=1)))
plot(gg)

## ----plot_IGKC-----------------------------------------------------------
PBMC_DR$IGKC = sc$toc['IGKC',rownames(PBMC_DR)]
gg = ggplot(PBMC_DR,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)

## ----sanity_check--------------------------------------------------------
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
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes),clusters=FALSE)

## ----visNullMatrix-------------------------------------------------------
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)

## ----calcNullMatrixWithClustering----------------------------------------
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes))
plotMarkerMap(sc,geneSet=igGenes,DR=PBMC_DR,useToEst=useToEst)

## ----calcContamination---------------------------------------------------
sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst)

## ----viewCont------------------------------------------------------------
head(sc$metaData)

## ----decontaminate-------------------------------------------------------
out = adjustCounts(sc)

## ----mostZeroed----------------------------------------------------------
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

## ----writeOut,eval=FALSE-------------------------------------------------
#  DropletUtils:::write10xCounts('./strainedCounts',out)

