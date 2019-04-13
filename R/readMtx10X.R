#' Load the 10X data matrix
#'
#' Loads 10X cellranger output from multiple directories and merge.
#'
#' @param Vector of directories containing matrix.mtx(.gz) files to load.
#' @param chNames Names of each experiment to prepend to barcodes to make them unique.
#' @param isV3 Are these cellranger v3.0.0 or greater outputs?  Can be a vector.
#' @return Matrix containing all data.
#' @importFrom Matrix readMM
readMtx10X = function(dataDirs,chNames=names(dataDirs),isV3=FALSE){
  #Ensure unique names
  if(is.null(chNames))
    chNames = paste0('Channel',seq_along(dataDirs))
  if(length(isV3)==1)
    isV3 = rep(isV3,length(dataDirs))
  out = list()
  for(i in seq_along(dataDirs)){
    #Load the raw matrix
    mtx = readMM(file.path(dataDirs[i],ifelse(isV3[i],'matrix.mtx.gz','matrix.mtx')))
    #Get row names.
    gns = read.table(file.path(dataDirs[i],ifelse(isV3[i],'features.tsv.gz','genes.tsv')),sep='\t',header=FALSE)
    rownames(mtx) = paste0(gns[,1],'___',gns[,2])
    #Get column names.
    noms = read.table(file.path(dataDirs[i],ifelse(isV3[i],'barcodes.tsv.gz','barcodes.tsv')),sep='\t',header=FALSE)
    colnames(mtx) = paste0(chNames[i],'___',noms[,1])
    out[[chNames[i]]] = mtx
  }
  do.call(cbind,out)
}
