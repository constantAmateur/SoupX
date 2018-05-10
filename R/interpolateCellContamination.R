#' Get cell specific rho from binned values
#'
#' Interpolates rho at each cell based on each cell's total number of UMIs from estimates in bins.
#'
#' @export
#' @param sc A \code{SoupChannel} or \code{SoupChannelList} object.
#' @param channelName The name of the channel to use if \code{sc} is a \code{SoupChannelList} object.
#' @param useGlobal If TRUE, set all cells to the global estimate of the contamination fraction.
#' @return A modified version of \code{sc} with cell estimates of contamination stored in \code{rhos}.  If \code{sc} is a \code{SoupChannelList}, the named channel is updated in this way.
interpolateCellContamination = function(sc,channelName,useGlobal=FALSE){
  if(is(sc,'SoupChannelList')){
    sc$channels[[channelName]] = interpolateCellContamination(sc$channels[[channelName]],useGlobal=useGlobal)
    return(sc)
  }
  if(!is(sc,'SoupChannel'))
    stop("scl must be an object of type SoupChannel or SoupChannelList")
  if(is.null(sc$rhoGrouped))
    stop("Run calculateContaminationFraction first.")
  #If we have limited values, just use the global one
  if(useGlobal | nrow(sc$rhoGrouped)<4){
    sc$rhos = (rep(sc$rhoGrouped['Global','est'],length(sc$nUMIs)))
  }else{
    #Exclude the global estimate
    rhos = sc$rhoGrouped[rownames(sc$rhoGrouped)!='Global',]
    fit = approx(rhos$nUMIs,rhos$est,nUMIs,rule=2)
    sc$rhos = fit$y
  }
  #Add names
  names(sc$rhos) = colnames(sc$toc)
  sc
}
