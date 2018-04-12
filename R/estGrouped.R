#' Estimates contamination from a group of cells
#'
#' Intended for internal use only.  Given expected and observed soup counts for a group of cells and a vector defining how to group, returns the maximum likelihood estimator of the contamination fraction for each group.
#'
#' @param obsCnts The observed number of counts that can confidently be assigned to background for each cell.
#' @param expCnts The expected number of counts for each cell if the cell were to contain nothing but background.
#' @param nUMIs The total UMIs for each cell.
#' @param grpFac A vector defining how cells are to be grouped together for calculating rho.
#' @param logRho Should a log10 transformation be applied to the estimate of rho?
#' @return A data.frame with estimates of rho by group, along with confidence bounds and some meta-data.
estGrouped = function(obsCnts,expCnts,nUMIs,grpFac,logRho=FALSE){
  #Drop any 0-groups from factor
  grpFac = factor(grpFac,levels=unique(as.character(grpFac)))
  #Calculate the quantifications
  oCnts = sapply(split(obsCnts,grpFac),sum)
  grpCnts = sapply(split(nUMIs,grpFac),sum)
  oCnts = estRateLims(oCnts,grpCnts)*grpCnts
  eCnts = sapply(split(expCnts,grpFac),sum)
  #Estimate and format for output
  df = oCnts/eCnts
  if(logRho)
    df = log10(df)
  #Record observed and expected counts
  df$obsSoupCnts = oCnts$est
  df$expSoupCnts = eCnts
  df$nUMIs = as.numeric(grpCnts/table(grpFac)[names(grpCnts)])
  df$isLogged=logRho
  return(df)
}
