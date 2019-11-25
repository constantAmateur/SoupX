#' Load a collection of 10X data-sets
#'
#' Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
#'
#' @export
#' @param dataDir Top level cellranger output directory (the directory that contains the "raw_gene_bc_matrices" folder).
#' @param cellIDs Barcodes of droplets that contain cells.  If NULL, use the default cellranger set.
#' @param channelName The name of the channel to store.  If NULL set to either \code{names(dataDir)} or \code{dataDir} is no name is set.
#' @param readArgs A list of extra parameters passed to \code{Seurat::Read10X}.
#' @param includeFeatures If multiple feature types are present, keep only the types mentioned here and collapse to a single matrix.
#' @param verbose Be verbose?
#' @param ... Extra parameters passed to \code{SoupChannel} construction function.
#' @return A SoupChannel object containing the count tables for the 10X dataset.
#' @seealso SoupChannel estimateSoup
#' @importFrom Seurat Read10X
#' @importFrom utils read.csv
load10X = function(dataDir,cellIDs=NULL,channelName=NULL,readArgs=list(),includeFeatures=c('Gene Expression'),verbose=TRUE,...){
  #Work out which version we're dealing with
  isV3 = dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))
  tgt = file.path(dataDir,
                  ifelse(isV3,'raw_feature_bc_matrix','raw_gene_bc_matrices'))
  #Add the reference genome for the non-V3 ones
  if(!isV3)
    tgt = file.path(tgt,list.files(tgt))
  if(verbose)
    message(sprintf("Loading raw count data"))
  dat = do.call(Read10X,c(list(data.dir=tgt),readArgs))
  if(verbose)
    message(sprintf("Loading cell-only count data"))
  if(!is.null(cellIDs)){
    #Do the same sripping that Seurat does on IDs
    if(all(grepl('\\-1$',cellIDs)))
      cellIDs = gsub('\\-1$','',cellIDs)
    #Check we have the IDs
    if(!all(cellIDs %in% colnames(dat)))
      stop("Not all supplied cellIDs found in raw data.")
    datCells = dat[,match(cellIDs,colnames(dat))]
  }else{
    #Work out which ones contain cells
    tgt = file.path(dataDir,
                    ifelse(isV3,'filtered_feature_bc_matrix','filtered_gene_bc_matrices'))
    if(!isV3)
      tgt = file.path(tgt,list.files(tgt))
    datCells = do.call(Read10X,c(list(data.dir=tgt),readArgs))
    #If it's a list of multiple types, have to decide what to include and collapse to one matrix.
    if(is.list(dat)){
      dat = do.call(rbind,dat[includeFeatures])
      datCells = do.call(rbind,datCells[includeFeatures])
    }
  }
  if(verbose)
    message(sprintf("Loading extra analysis data where available"))
  #Get the cluster annotation if available
  mDat = NULL
  tgt = file.path(dataDir,'analysis','clustering','graphclust','clusters.csv')
  if(file.exists(tgt)){
    clusters = read.csv(tgt)
    mDat = data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
  }
  #Add fine grained clusters too if present
  tgt = file.path(dataDir,'analysis','clustering','kmeans_10_clusters','clusters.csv')
  if(file.exists(tgt)){
    clusters = read.csv(tgt)
    mDat$clustersFine = clusters$Cluster
  }
  #Get tSNE if available and point to it
  tgt = file.path(dataDir,'analysis','tsne','2_components','projection.csv')
  if(file.exists(tgt)){
    tsne = read.csv(tgt)
    if(is.null(mDat)){
      mDat = data.frame(tSNE1=tsne$TSNE.1,tSNE2=tsne$TSNE.2,row.names=tsne$Barcode)
    }else{
      mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat),tsne$Barcode)]
      mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat),tsne$Barcode)]
    }
    DR = c('tSNE1','tSNE2')
  }else{
    DR=NULL
  }
  #Ensure rownames of metadata match column names of counts
  if(!is.null(mDat) && any(rownames(mDat)!=colnames(datCells))){
    rownames(mDat) = gsub('-1$','',rownames(mDat))
    if(any(rownames(mDat)!=colnames(datCells)))
      stop("Error matching meta-data to cell names.")
  }
  #Get a name for the channel
  if(is.null(channelName))
    channelName = ifelse(is.null(names(dataDir)),dataDir,names(dataDir))
  channel = SoupChannel(tod = dat,
                        toc = datCells,
                        metaData = mDat,
                        channelName = channelName,
                        dataDir = dataDir,
                        dataType='10X',
                        isV3=isV3,
                        DR=DR,
                        ...
                        )
  return(channel)
}
