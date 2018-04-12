#' Get cell specific rho from binned values
#'
#' Interpolates rho at each cell based on each cell's total number of UMIs from estimates in bins.
#'
#' @export
#' @param rhos The output of \code{\link{calculateContaminationFraction}} for this set of cells.
#' @param nUMIs The total number of UMIs for each cell.
#' @return A vector of the same length as nUMIs with a cell-specific estimate of rho.
interpolateCellContamination = function(rhos,nUMIs){
  #If we have limited values, just use the global one
  if(nrow(rhos)<4)
    return(rep(rhos['Global','est'],length(nUMIs)))
  #Exclude the global estimate
  rhos = rhos[rownames(rhos)!='Global',]
  fit = approx(rhos$nUMIs,rhos$est,nUMIs,rule=2)
  return(fit$y)
}
