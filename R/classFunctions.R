#' Construct a SoupChannel object
#'
#' Creates a SoupChannel object that contains everything related to the soup estimation of a single channel.
#'
#' @export
#' @param tod Table of droplets.  A matrix with columns being each droplet and rows each gene.
#' @param toc Table of counts.  Just those columns of \code{tod} that contain cells.
#' @param metaData Meta data pertaining to the cells.  Optional.  Must be a data-frame with rownames equal to column names of \code{toc}.
#' @param calcSoupProfile By default, the soup profile is calculated using \code{\link{estimateSoup}} with default values.  If you want to do something other than the defaults, set this to \code{FALSE} and call \code{\link{estimateSoup}} manually.
#' @param ... Any other named parameters to store.
#' @return A SoupChannel object.
#' @importFrom Matrix colSums
#' @examples
#' \dontrun{
#' #Default calculates soup profile
#' sc = SoupChannel(tod,toc)
#' names(sc)
#' #This can be suppressed
#' sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
#' names(sc)
#' #And the soup profile calculated with non-default values
#' sc = estimateSoup(sc,soupRange=c(10,50))
#' #Or added manually
#' soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),counts=rowSums(toc))
#' sc = setSoupProfile(sc,soupProf)
#' }
#' @seealso SoupChannelList estimateSoup
SoupChannel = function(tod,toc,metaData=NULL,calcSoupProfile=TRUE,...){
  if(!is.null(metaData) & !all(sort(colnames(toc))==sort(rownames(metaData))))
    stop("Rownames of metaData must match column names of table of counts.")
  #Munge everything into a list
  out = list(tod=tod,toc=toc)
  out = c(out,list(...))
  #Create the metadata object
  out$metaData = data.frame(row.names=colnames(toc),
                            nUMIs = colSums(toc)
                            )
  #Merge in supplied metaData if it's present
  if(!is.null(metaData)){
    #Drop nUMIs if it exists
    metaData = metaData[,colnames(metaData)!='nUMIs',drop=FALSE]
    out$metaData = cbind(out$metaData,metaData)
  }
  #Get the droplet UMIs as well, as that's a useful thing to have
  out$nDropUMIs = colSums(tod)
  class(out) = c('list','SoupChannel')
  #Estimate the soup
  if(calcSoupProfile)
    out = estimateSoup(out)
  return(out)
}

#' Print method for SoupChannel
#'
#' Prints a summary of a SoupChannel object.
#' 
#' @export
#' @param x A SoupChannel object.
#' @param ... Currently unused.
print.SoupChannel = function(x,...) {
  message(sprintf("Channel with %d genes and %d cells",nrow(x$toc),ncol(x$toc)))
}

