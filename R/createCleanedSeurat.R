#' Convert multiple strained channels into a Seurat object
#'
#' Constructs a Seurat object from the channels given and normalises the result.
#'
#' @export
#' @param cleanedChannels A list of channels that have been run through \code{\link{strainChannel}}.
#' @param channelPrefixes A vector of the same length as \code{cleanedChannels} giving a prefix to add to each cell before merging.
#' @param ... Extra parameters passed to \code{CreateSeuratObject}.
#' @return A normalised Seurat object.
createCleanedSeurat = function(cleanedChannels,scale.factor=1e4,channelCellPrefix='',...){
  #Fix prefix object if necessary
  if(length(channelCellPrefix)==1)
    channelCellPrefix = rep(channelCellPrefix,length(cleanedChannels))
  #Create big TOC with raw data to initialise
  bigTOC = do.call(cbind,lapply(seq_along(cleanedChannels),function(e){
                                tmp = cleanedChannels[[e]]$rawCounts
                                colnames(tmp) = paste0(channelCellPrefix[e],colnames(tmp))
                                tmp}))
  srat = CreateSeuratObject(bigTOC,...)
  #Create the cleaned replacement table
  cleaned = do.call(cbind,lapply(seq_along(cleanedChannels),function(e){
                                tmp = cleanedChannels[[e]]$correctedProfile
                                colnames(tmp) = paste0(channelCellPrefix[e],colnames(tmp))
                                tmp}))
  #Do the normalisation
  cleaned = log(1+cleaned*scale.factor)
  #Now cram into the normalisation slot.  This is ugggggggly
  params = list(object=srat,assay.type='RNA',normalization.method='LogNormalize',scale.factor=scale.factor,display.progress=TRUE)
  srat = Seurat:::SetCalcParams(object = srat,calculation='NormalizeData',...=params)
  srat = SetAssayData(object=srat,assay.type='RNA',slot='data',new.data=cleaned)
  #Update the meta-data to account for correction
  srat@meta.data$nGene = colSums(srat@data>0)
  return(srat)
}


