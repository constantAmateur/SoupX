#' Get expression profile of soup
#'
#' Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.
#'
#' @export
#' @param tod Table of UMIs for all droplets.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @return A data.frame with the soup profile and confidence limits for all genes.
estimateSoup = function(tod,soupRange=c(0,10)){
  nUMIs = colSums(tod)
  w = which(nUMIs>soupRange[1] & nUMIs<soupRange[2])
  estRateLims(rowSums(tod[,w,drop=FALSE]),sum(tod[,w]))
}
