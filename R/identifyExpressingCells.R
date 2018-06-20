#' Identify expressing cells 
#'
#' This function identifies which cells truly express different sets of genes and which only appear to express it due to ambient mRNA contamination.  A list containing genes that can be expected to be non-expressed in at least one cell must be given.  That is, this list contains sets of genes (e.g. Haemoglobin genes) that can be assumed to be completely absent from certain cells.  For example, haemoglobin genes should not be expressed in anything that is not a red blood cell.
#'
#'  The different methods to determine which cells to exclude which genes for are:
#'   pCut - Exclude any cell where we can reject the null hypothesis (Poisson test) that the contamination fraction is less than or equal to 1 for a set of genes at the \code{exCut} level.
#'   thresh - Exclude any cell for which the estimate for contamination fraction exceeds \code{exCut}.
#'
#' @export 
#' @param sc A SoupChannel or SoupChannelList object.  If a SoupChannelList object, channelName must be given and must be the name of a channel.
#' @param channelName The name of a channel to use if \code{sc} is a SoupChannelList object.
#' @param nonExpressedGeneList A list containing sets of genes which can be assumed to be non-expressed in a subset of cells (see details).
#' @param excludeMethod Which method to use to exclude cells (see details).
#' @param exCut Cut-off used by \code{excludeMethod} (see details).
#' @return A boolean matrix with cells that definitively express the set of genes given by the row name set to FALSE.  Usually passed to \code{\link{calculateContaminationFraction}} via the \code{useToEst} argument.
#' @seealso \code{\link{calculateContaminationFraction}} 
identifyExpressingCells = function(sc,channelName,nonExpressedGeneList,excludeMethod=c('pCut','thresh'),exCut=0.05) {
  if(is(sc,'SoupChannelList')){
    if(!(channelName %in% names(sc$channels)))
      stop("sc is a SoupChannelList object, but channelName is not the name of a channel within this object")
    return(identifyExpressingCells(sc$channels[[channelName]],'',nonExpressedGeneList=nonExpressedGeneList,excludeMethod=excludeMethod,exCut=exCut))
  }
  excludeMethod = match.arg(excludeMethod)
  #Collapse genes in each marker list down into one measure
  sFracs = do.call(rbind,lapply(nonExpressedGeneList,function(e) estRateLims(sum(sc$soupProfile[e,'cnts']),sum(sc$soupProfile$cnts))))
  #Do the same with cells, we don't care about other genes for estimation (other than the total counts)
  dat = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(sc$toc[e,,drop=FALSE])))
  if(is.null(dim(dat)))
    dat = matrix(dat,nrow=1,dimnames=list(names(nonExpressedGeneList),names(dat)))
  #############################
  # Create cell mask for soup
  #For each marker set, work out which cells to exclude
  useToEst = t(sapply(names(nonExpressedGeneList),function(e) {
                  if(excludeMethod=='thresh'){
                    return((dat[e,]/(sc$nUMIs*sFracs[e,'est']))<exCut)
                  }else{
                    qVals = ppois(dat[e,]-1,sc$nUMIs*sFracs[e,'est'],lower.tail=FALSE)
                    qVals = p.adjust(qVals,method='BH')
                    ifelse(qVals<exCut,FALSE,TRUE)
                  }
            }))
  return(useToEst)
}
 
