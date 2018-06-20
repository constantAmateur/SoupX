#' Estimates which droplets contain a cell
#'
#' This function attempts to identify which droplets contain cells based on the level of contamination in each cell.  For this to work well, there must be sufficient power to estimate the contamination fraction in each droplet, which is usually not the case unless spike-ins have been used.  To determine which droplets contain cells, a binomial test (with FDR correction) is performed against the NULL hypothesis that the contamination rate is equal to or greater than \code{contCut}.
#'
#' @export
#' @param sc A SoupChannel object for which the table of droplets has been retained.  That is, \code{\link{estimateSoup}} was run with \code{keepDroplets=TRUE}.
#' @param nonExpressedGeneList Genes to use for estimating contamination fraction (see \code{\link{calculateContaminationFraction}})
#' @param useToEst Which cells to use which sets of genes to estimate contamination (see \code{\link{calculateContaminationFraction}})
#' @param minUMIs Don't test any droplet with fewer than this many UMIs.
#' @param minGenes Don't test any droplet with fewer than this many detected genes.
#' @param contCut The maximum allowable contamination fraction for a droplet containing a cell.
#' @return A data.frame with contamination fraction estimates for each droplet and a p-value and FDR for rejecting the NULL hypothesis and accepting this droplet as a cell.
findCells = function(sc,nonExpressedGeneList,useToEst=NULL,minUMIs=100,minGenes=50,contCut=0.3) {
  nUMIs = colSums(sc$tod)
  nGenes = colSums(sc$tod>0)
  w = which(nUMIs>minUMIs & nGenes > minGenes)
  #Estimate rho at the individual cell level
  sc$toc = sc$tod[,w]
  rhos = calculateContaminationFraction(sc,nonExpressedGeneList,useToEst=useToEst,cellGroups=colnames(sc$tod)[w])$rhoGrouped
  #Drop global
  rhos = rhos[rownames(rhos)!='Global',]
  #Test the NULL hypothesis is that the contamination fraction is greater than or equal to x
  pVals = ppois(rhos$obsSoupCnts,rhos$expSoupCnts*contCut,lower.tail=TRUE)
  qVals = p.adjust(pVals,method='BH')
  rhos$pVals = pVals
  rhos$FDR=qVals
  return(rhos)
}
