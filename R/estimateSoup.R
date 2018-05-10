#' Get expression profile of soup
#'
#' Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.  If a SoupChannelList is provided, the soup profile is estimated for each channel and stored.
#'
#' @export
#' @param scl A \code{SoupChannelList}.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets Storing the full table of counts for all droplets uses a lot of space and is really only used to estimate the soup profile.  Therefore, it is dropped after the soup profile has been estimated unless this is set to \code{TRUE}.
#' @return A modified version of \code{scl} with an extra \code{soupProfile} entry containing a data.frame with the soup profile and confidence limits for all genes.
estimateSoup = function(scl,soupRange=c(0,10),keepDroplets=FALSE){
  if(!is(scl,'SoupChannelList'))
    stop("scl must be a SoupChannelList object.")
  #Estimate the soup for each channel
  for(i in seq_along(scl$channels)){
    nUMIs = colSums(scl$channels[[i]]$tod)
    w = which(nUMIs > soupRange[1] & nUMIs < soupRange[2])
    scl$channels[[i]]$soupProfile = estRateLims(rowSums(scl$channels[[i]]$tod[,w,drop=FALSE]),sum(scl$channels[[i]]$tod[,w]))
    if(!keepDroplets)
      scl$channels[[i]]$tod=NULL
  }
  #Create a convenient table to look up soup expression
  scl$soupMatrix = do.call(cbind,lapply(scl$channels,function(e) e$soupProfile[,'est']))
  rownames(scl$soupMatrix) = rownames(scl$toc)
  return(scl)
}
