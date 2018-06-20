#' Convert multiple strained channels into a Seurat object
#'
#' Constructs a Seurat object from the channels given and normalises the result.
#'
#' @export
#' @param scl A \code{SoupChannelList} object on which \code{\link{strainCells}} has been run.
#' @param scale.factor The scale factor used in normalising data to log(1+scale.factor*x)
#' @param ... Extra parameters passed to \code{CreateSeuratObject}.
#' @return A normalised Seurat object.
#' @importFrom Seurat CreateSeuratObject SetAssayData
createCleanedSeurat = function(scl,scale.factor=1e4,...){
  if(!is(scl,'SoupChannelList'))
    stop("scl must be a SoupChannelList object")
  #Initialise with table of counts
  srat = CreateSeuratObject(scl$toc,...)
  #Create the cleaned replacement table
  cleaned = log(1+scl$strainedExp*scale.factor)
  #Now cram into the normalisation slot.  This is ugggggggly
  params = list(object=srat,assay.type='RNA',normalization.method='LogNormalize',scale.factor=scale.factor,display.progress=TRUE)
  #Pulled from Seurat:::SetCalcParams
  srat@calc.params[['NormalizeData']] = params
  srat@calc.params[['NormalizeData']]$object = NULL
  srat@calc.params[['NormalizeData']]$object2 = NULL
  srat@calc.params[['NormalizeData']]$time = Sys.time()
  #srat = Seurat:::SetCalcParams(object = srat,calculation='NormalizeData',...=params)
  srat = SetAssayData(object=srat,assay.type='RNA',slot='data',new.data=cleaned)
  #Update the meta-data to account for correction
  srat@meta.data$nGene = colSums(srat@data>0)
  return(srat)
}
