#' Estimates background contamination from a single sample
#'
#' Estimates the background contamination profile, cell-specific contamination fraction and the expression profile for each cell without contamination.
#'
#' @export
#' @param tod Table of counts (UMIs) where columns are samples and rows are genes.  Must include all droplets, not just those containing cells (e.g., raw_gene_bc_matricies for 10X data)
#' @param cellIdxs Indices for the columns of \code{toc} that contain cells.
#' @param nonExpressedGeneList Which genes to use to estimate soup (see \code{\link{calculateContaminationFraction}})
#' @param soupRange Which droplets to estimate the soup from (see \code{\link{estimateSoup}})
#' @param ... Extra arguments passed to \code{\link{calculateContaminationFraction}}.
#' @importFrom Matrix colSums
strainChannel = function(tod,cellIdxs,nonExpressedGeneList,soupRange=c(0,10),...){
  #Estimate the soup
  soupProfile = estimateSoup(tod,soupRange)
  #Calculate rho for each cell
  rhos = calculateContaminationFraction(tod[,cellIdxs],soupProfile,nonExpressedGeneList,...)
  #Interpolate rho to each cell
  cellRhos = interpolateCellContamination(rhos,colSums(tod[,cellIdxs]))
  #Now de-contaminate and return estimate
  trueCellExpression = strainCells(tod[,cellIdxs],cellRhos,soupProfile)
  #And also calculate the ratio of the observed counts to the soup.  Kept un-logged so it remains sparse
  #We should really convert a bunch of 0s to NaNs, but don't do this to save space.
  expressionRatio = t(t(toc)/colSums(toc))
  expressionRatio@x = expressionRatio@x/soupProfile[expressionRatio@i+1,'est']
  #Now return everything that we've calculated
  return(list(contaminationFractions = rhos, 
              imputedCellContaminationFractions = cellRhos, 
              rawCounts = tod[,cellIdxs],
              correctedProfile = trueCellExpression, 
              ratioMatrix = expressionRatio,
              soupProfile = soupProfile,
              nonExpressedGeneList = nonExpressedGeneList))
}


