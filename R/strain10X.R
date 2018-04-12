#' Estimates background contamination and removes it from a 10X channel
#'
#' Takes the definition of cells from the "filtered counts" folder, then runs the standard workflow to estimate and remove contamination.
#'
#' @export
#' @param dataDirs Path to the 10X Cell Ranger output.  Should be the top level folder containing the html report.  Can be a vector, in which case each is processed in turn with the all other parameters unchanged.
#' @param nonExpressedGeneList Which genes to use to estimate soup (see \code{\link{calculateContaminationFraction}})
#' @param soupRange Which droplets to estimate the soup from (see \code{\link{estimateSoup}})
#' @param ... Extra arguments passed to \code{\link{calculateContaminationFraction}}.
#' @importFrom Seurat Read10X
#' @return The same as \code{\link{strainChannel}}
strain10X = function(dataDirs,nonExpressedGeneList,soupRange=c(0,10),...){
  channels=list()
  for(dataDir in dataDirs){
    #Get reference
    ref = list.files(file.path(dataDir,'raw_gene_bc_matrices'))
    #Load the 10X data
    tod = Read10X(file.path(dataDir,'raw_gene_bc_matrices',ref))
    #Get the barcodes that cell ranger considers to contain cells
    cells = read.delim(file.path(dataDir,'filtered_gene_bc_matrices',ref,'barcodes.tsv'),sep='\t',header=FALSE)
    cells = gsub('-1','',cells[,1])
    #Get the index in the big table
    cellIdxs = match(cells,colnames(tod))
    #Now process this channel in the usual way
    channels[[dataDir]] = strainChannel(tod,cellIdxs,nonExpressedGeneList,soupRange=soupRange,...)
    #Record the reference genome
    channels[[dataDir]]$refGenome = ref
  }
  if(length(channels)==1)
    return(channels[[1]])
  return(channels)
}


