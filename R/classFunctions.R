#' Construct a SoupChannel object
#'
#' Creates a SoupChannel object that contains everything related to the soup estimation of a single channel.
#'
#' @export
#' @param tod Table of droplets.  A matrix with columns being each droplet and rows each gene.
#' @param toc Table of counts.  Just those columns of \code{tod} that contain cells.
#' @param channelName A name for this channel.
#' @param ... Any other named parameters to store.
#' @return A SoupChannel object.
SoupChannel = function(tod,toc,channelName,...){
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
  class(out) = c('list','SoupChannel')
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
SoupChannelList = function(channels,...){
  if(any(duplicated(sapply(channels,function(e) e$channelName))))
    stop("Duplicate channel names found.  Please give each channel a unique name before continuing.")
  #Ensure consistent naming
  names(channels) = sapply(channels,function(e) e$channelName)
  scl = list(channels=channels)
  scl$toc = do.call(cbind,lapply(channels,function(e) e$toc))
  scl$nUMIs = colSums(scl$toc)
  scl = c(scl,list(...))
  class(scl) = c('list','SoupChannelList')
  scl
}


#' Print method for SoupChannel
#'
#' Prints a summary of a SoupChannel object.
#' 
#' @export
#' @param sc A SoupChannel object.
print.SoupChannel = function(sc) {
  message(sprintf("Channel named %s with %d genes and %d cells",sc$channelName,nrow(sc$toc),ncol(sc$toc)))
}

#' Print method for SoupChannelList
#'
#' Prints a summary of a SoupChannelList object.
#' 
#' @export
#' @param scl A SoupChannelList object.
print.SoupChannelList = function(scl) {
  message(sprintf("%d Channels with %d genes and %d cells in total.",length(scl$channels),nrow(scl$toc),ncol(scl$toc)))
}
