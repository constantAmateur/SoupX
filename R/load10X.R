#' Load a collection of 10X data-sets
#'
#' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'
#' @export
#' @param dataDirs Vector of top level cellranger output directories (the directory that contains the "raw_gene_bc_matrices" folder).
#' @param channelNames To make droplet barcodes unique across experiment, each channel needs its own unique label.  If not given, this is set numerically.
#' @param dropEns By default, genes are named in the format <ENSEMBL_ID>___<Gene_Symbol>.  If dropEns is set to TRUE, the ENESEMBL ID is dropped from the name.
#' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#' @return A SoupChannelList object containing the count tables for each 10X dataset.
#' @seealso SoupChannel SoupChannelList estimateSoup
load10X = function(dataDirs,channelNames=names(dataDirs),dropEns=TRUE,...){
  if(is.null(channelNames))
    channelNames = sprintf('Channel%d',seq_along(dataDirs))
  channels = list()
  for(i in seq_along(dataDirs)){
    message(sprintf("Loading data for 10X channel %s from %s",channelNames[i],dataDirs[i]))
    dataDir = dataDirs[i]
    #Work out which version we're dealing with
    ref = file.path(dataDir,'raw_gene_bc_matrices')
    isV3=FALSE
    if(!dir.exists(ref)){
      if(dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))){
        isV3=TRUE
      }else{
        stop(sprintf("Invalid cellranger directory %s",dataDir))
      }
    }
    #Load matrix and cell IDs
    if(isV3){
      ref=NA
      tod = readMtx10X(file.path(dataDir,'raw_feature_bc_matrix'),channelNames[i],isV3=TRUE)
      cells = read.delim(file.path(dataDir,'filtered_feature_bc_matrix','barcodes.tsv.gz'),sep='\t',header=FALSE)
      cells = paste0(channelNames[i],'___',cells[,1])
    }else{
      ref = list.files(ref)
      #Load the 10X data
      tod = readMtx10X(file.path(dataDir,'raw_gene_bc_matrices',ref),channelNames[i],isV3=FALSE)
      #Get the barcodes that cell ranger considers to contain cells
      cells = read.delim(file.path(dataDir,'filtered_gene_bc_matrices',ref,'barcodes.tsv'),sep='\t',header=FALSE)
      cells = paste0(channelNames[i],'___',cells[,1])
    }
    if(dropEns)
      rownames(tod) = make.unique(gsub('.*___','',rownames(tod)))
    #Get the index in the big table
    cellIdxs = match(cells,colnames(tod))
    #Get the cluster annotation if available
    mDat = NULL
    tgt = file.path(dataDir,'analysis','clustering','graphclust','clusters.csv')
    if(file.exists(tgt)){
      clusters = read.csv(tgt)
      mDat = data.frame(clusters=clusters$Cluster,row.names=paste0(channelNames[i],'___',clusters$Barcode))
    }
    #Add fine grained clusters too if present
    tgt = file.path(dataDir,'analysis','clustering','kmeans_10_clusters','clusters.csv')
    if(file.exists(tgt)){
      clusters = read.csv(tgt)
      mDat$clustersFine = clusters$Cluster
    }
    #Get tSNE if available
    tgt = file.path(dataDir,'analysis','tsne','2_components','projection.csv')
    if(file.exists(tgt)){
      tsne = read.csv(tgt)
      cellIDs = paste0(channelNames[i],'___',tsne$Barcode)
      if(is.null(mDat)){
        mDat = data.frame(tSNE1=tsne$TSNE.1,tSNE2=tsne$TSNE.2,row.names=cellIDs)
      }else{
        mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat),cellIDs)]
        mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat),cellIDs)]
      }
    }
    channels[[channelNames[i]]] = SoupChannel(tod,tod[,cellIdxs,drop=FALSE],channelName=channelNames[i],metaData=mDat,ref=ref,path=dataDir,dataType='10X',isV3=isV3,...)
  }
  channels = SoupChannelList(channels)
  return(channels)
}
