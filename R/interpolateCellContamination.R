#' Get cell specific rho from binned values
#'
#' Interpolates rho at each cell based on each cell's total number of UMIs from estimates in bins.  By default, uses a lowess fit.  But options for linear interpolation or simply using the global value are also available.
#' 
#'If it is not very well known which cells truly express the genes used to estimate the contamination (e.g. which cells are red blood cells when using haemoglobin to estimate contamination), the estimate of rho is often "spikey".  That is, there will be clear overall trend, but bins that erroneously contain cells that should not have been used to estimate rho but have will have a very high contamination fraction.  In such cases it is better to interpolate using either lowess or the global estimate.
#'
#'
#' @export
#' @param sc A \code{SoupChannel} or \code{SoupChannelList} object.
#' @param channelName The name of the channel to use if \code{sc} is a \code{SoupChannelList} object.
#' @param interpolationMethod Method to use to interpolate.  \code{lowess} uses a constrained monotonically decreasing lowess fit, \code{lowess} uses a lowess fit without constraints, \code{linear} uses linear interpolation and \code{fixed} set all cells to have the global contamination value, which is set by \code{fixedContaminationValue}.
#' @param fixedContaminationValue The fixed contamination value to used when \code{interpolationMethod} is set to "fixed".  If set to NA, use the global contamination estimate.
#' @param useGlobal Sets interpolationMethod to \code{global}.  Retained for compatibility, use \code{interpolationMethod} directly instead.
#' @param ... Parameters to pass to interpolation method.
#' @return A modified version of \code{sc} with cell estimates of contamination stored in \code{rhos}.  If \code{sc} is a \code{SoupChannelList}, the named channel is updated in this way.
interpolateCellContamination = function(sc,channelName,interpolationMethod = c('decreasingLowess','lowess','linear','fixed'),fixedContaminationValue=NA,useGlobal=FALSE,...){
  interpolationMethod = match.arg(interpolationMethod)
  if(is(sc,'SoupChannelList')){
    sc$channels[[channelName]] = interpolateCellContamination(sc$channels[[channelName]],interpolationMethod=interpolationMethod,fixedContaminationValue=fixedContaminationValue,useGlobal=useGlobal)
    return(sc)
  }
  if(useGlobal){
    interpolationMethod='fixed'
    fixedContaminationValue=NA
  }
  if(!is(sc,'SoupChannel'))
    stop("scl must be an object of type SoupChannel or SoupChannelList")
  if(is.null(sc$rhoGrouped))
    stop("Run calculateContaminationFraction first.")
  #If we have limited values, just use the global one
  if(interpolationMethod=='fixed' | nrow(sc$rhoGrouped)<4){
    #Get or set the fixed rho
    if(is.na(fixedContaminationValue))
      fixedContaminationValue = sc$rhoGrouped['Global','est']
    sc$rhos = (rep(fixedContaminationValue,length(sc$nUMIs)))
  }else if(interpolationMethod=='linear'){
    #Exclude the global estimate
    rhos = sc$rhoGrouped[rownames(sc$rhoGrouped)!='Global',]
    fit = approx(rhos$nUMIs,rhos$est,sc$nUMIs,rule=2)
    sc$rhos = fit$y
  }else if(interpolationMethod=='lowess'){
    #Use lowess fit
    w = rownames(sc$rhoGrouped)!='Global' & sc$rhoGrouped$nUMIs>0 & is.finite(sc$rhoGrouped$est) & !is.na(sc$rhoGrouped$est)
    rhos = sc$rhoGrouped[w,]
    fit = with(rhos,lowess(est ~ nUMIs,...))
    #Now need to interpolate on the smoothed values
    fitCells = approx(fit$x,fit$y,sc$nUMIs,rule=2)
    sc$rhos = fitCells$y
  }else{
    #Use monotonic decreasing lowess
    w = rownames(sc$rhoGrouped)!='Global' & sc$rhoGrouped$nUMIs>0 & is.finite(sc$rhoGrouped$est) & !is.na(sc$rhoGrouped$est)
    rhos = sc$rhoGrouped[w,]
    fit = with(rhos,lowess(est ~ nUMIs,...))
    ##Exclude anything that violates the monotonicity constraint
    #w = fit$y <= cummin(fit$y)
    ##And fill it in with a cubic spline.
    #fit$y[!w] = spline(fit$x[w],fit$y[w],xout=fit$x)$y[!w]
    #Anything that still doesn't fit, just interpolate over it.
    w = fit$y <= cummin(fit$y)
    if(sum(w)>1){
      fit$y[!w] = approx(fit$x[w],fit$y[w],fit$x,rule=2)$y[!w]
    }else{
      fit$y = rep(fit$y[w],length(fit$y))
    }
    #Finally do the linear interpolation
    if(length(fit$x)>1){
      fitCells = approx(fit$x,fit$y,sc$nUMIs,rule=2)
    }else{
      fitCells = rep(fit$y,sc$nUMIs)
    }
    sc$rhos = fitCells$y
  }
  #Make sure nothing is negative (can happen with monotonic interpolation)
  sc$rhos[sc$rhos<0] = 0
  #Add names
  names(sc$rhos) = colnames(sc$toc)
  sc
}
