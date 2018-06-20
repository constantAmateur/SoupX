#' Construct a SoupChannel object
#'
#' Creates a SoupChannel object that contains everything related to the soup estimation of a single channel.
#'
#' @export
#' @param tod Table of droplets.  A matrix with columns being each droplet and rows each gene.
#' @param toc Table of counts.  Just those columns of \code{tod} that contain cells.
#' @param channelName A name for this channel.
#' @param soupRange Which droplets to estimate soup from. Passed to \code{\link{estimateSoup}}.
#' @param keepDroplets Should we keep the full table of droplets?  Passed to \code{\link{estimateSoup}}.
#' @param ... Any other named parameters to store.
#' @return A SoupChannel object.
#' @seealso SoupChannelList estimateSoup
SoupChannel = function(tod,toc,channelName,soupRange=c(0,10),keepDroplets=FALSE,...){
  if(missing(channelName))
    channelName = 'UnknownChannel'
  #Make sure channelName is at the start of all channel names
  if(!all(gsub('___.*','',colnames(tod))==channelName))
    colnames(tod) = paste0(channelName,'___',colnames(tod))
  if(!all(gsub('___.*','',colnames(toc))==channelName))
    colnames(toc) = paste0(channelName,'___',colnames(toc))
  out = list(tod=tod,toc=toc,channelName=channelName)
  out = c(out,list(...))
  out$nUMIs = colSums(toc)
  #Get the droplet UMIs as well, as that's a useful thing to have
  out$nDropUMIs = colSums(tod)
  class(out) = c('list','SoupChannel')
  #Estimate the soup
  out = estimateSoup(out,soupRange=soupRange,keepDroplets=keepDroplets)
  out
}

#' Construct a SoupChannelList object
#'
#' Creates a SoupChannelList object that contains everything related to the soup estimation for a collection of channels, usually one experiment.
#'
#' @export
#' @param channels A list of \code{\link{SoupChannel}} objects.  Must be uniquely named.
#' @param ... Any other named parameter to store.
#' @return A SoupChannelList object.
#' @seealso SoupChannel estimateSoup
SoupChannelList = function(channels,...){
  if(any(duplicated(sapply(channels,function(e) e$channelName))))
    stop("Duplicate channel names found.  Please give each channel a unique name before continuing.")
  #Ensure consistent naming
  names(channels) = sapply(channels,function(e) e$channelName)
  scl = list(channels=channels)
  scl$toc = do.call(cbind,lapply(channels,function(e) e$toc))
  scl$nUMIs = do.call(c,lapply(channels,function(e) e$nUMIs))
  scl = c(scl,list(...))
  class(scl) = c('list','SoupChannelList')
  #Create summary of soup
  scl$soupMatrix = do.call(cbind,lapply(scl$channels,function(e) e$soupProfile[,'est']))
  rownames(scl$soupMatrix) = rownames(scl$toc)
  scl
}


#' Print method for SoupChannel
#'
#' Prints a summary of a SoupChannel object.
#' 
#' @export
#' @param x A SoupChannel object.
#' @param ... Currently unused.
print.SoupChannel = function(x,...) {
  message(sprintf("Channel named %s with %d genes and %d cells",x$channelName,nrow(x$toc),ncol(x$toc)))
}

#' Print method for SoupChannelList
#'
#' Prints a summary of a SoupChannelList object.
#' 
#' @export
#' @param x A SoupChannelList object.
#' @param ... Currently unused.
print.SoupChannelList = function(x,...) {
  message(sprintf("%d Channels with %d genes and %d cells in total.",length(x$channels),nrow(x$toc),ncol(x$toc)))
}
