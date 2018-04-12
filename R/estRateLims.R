#' Calculate limits on estimated binomial rate from data
#'
#' Internal use only.  Gets 95% confidence limits on an observation of x counts for something out of a total of n observations.
#'
#' @param x Number of positive observations.
#' @param n Out of this many.
#' @param conf Confidence interval
#' @param noms Row names to give to output.  If NULL, set to names(x).
#' @param sort Should we sort the output by est?
#' @return A data.frame giving estimated rate (x/n) and the 95% confidence limits.
estRateLims = function(x,n,conf=0.95,noms=NULL,sort=FALSE){
  if(is.null(noms))
    noms = names(x)
  #Coerce to integers
  x = as.integer(x)
  n = as.integer(n)
  alpha = (1-conf)/2
  p.L = ifelse(x==0,0,qbeta(alpha,x,n-x+1))
  p.U = ifelse(x==n,1,qbeta(1-alpha,x+1,n-x))
  x = (data.frame(est=x/n,lower=p.L,upper=p.U,cnts=x,total=n,row.names=noms))
  if(sort)
    return(x[order(x$est,decreasing=TRUE),])
  return(x)
}
