#' Get expression profile of soup
#'
#' This is usually called by \code{\link{SoupChannel}}, rather than directly by the user.  Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.
#'
#' @export
#' @param sc A \code{SoupChannel} object.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets Storing the full table of counts for all droplets uses a lot of space and is really only used to estimate the soup profile.  Therefore, it is dropped after the soup profile has been estimated unless this is set to \code{TRUE}.
#' @return A modified version of \code{sc} with an extra \code{soupProfile} entry containing a data.frame with the soup profile and confidence limits for all genes.
estimateSoup = function(sc,soupRange=c(0,100),keepDroplets=FALSE){
  if(!is(sc,'SoupChannel'))
    stop("sc must be a SoupChannel object.")
  #Estimate the soup 
  w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  sc$soupProfile = data.frame(row.names=rownames(sc$tod),
                              est = rowSums(sc$tod[,w,drop=FALSE])/sum(sc$tod[,w]),
                              counts = rowSums(sc$tod[,w,drop=FALSE]))
  #Saves a lot of space if we can drop the droplets now we're done with them
  if(!keepDroplets)
    sc$tod=NULL
  return(sc)
}
