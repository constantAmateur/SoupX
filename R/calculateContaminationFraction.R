#' Calculate the contamination fraction
#'
#' This function calculates the contamination fraction using the markers provided in groups of cells.  A list containing genes that can be expected to be non-expressed in at least one cell must be given.  That is, this list contains sets of genes (e.g. Haemoglobin genes) that can be assumed to be completely absent from certain cells.  For example, haemoglobin genes should not be expressed in anything that is not a red blood cell.  The contamination fraction is then calculated using these genes in groups of cells.
#'
#' As well as sets of non expressed genes, the user must indicate which cells these genes are non-expressed in.  Continuing with the haemoglobin gene example, this would mean the user must indicate which cells contain red blood cells.  This information is given by the parameter \code{useToEst} which should be a boolean matrix with columns representing cells in toc and rows representing groups of genes defined in \code{nonExpressedGeneList}.  If this is set to NULL, the code will attempt to estimate which cells truly express each set of genes in \code{nonExpressedGeneList} based on the method given in \code{cellExpressionMethod}.
#'
#' Unless a large list of non expressed genes can be given, there is usually insufficient power to estimate the contamination fraction separately for each cell.  Therefore, the user should supply a vector to \code{cellGroups} that can be used to group together cells when estimating the contamination fraction.  If set to NULL, the code will group together cells with a similar number of total UMIs in groups such that the expected number of soup counts in each group is roughly equal to \code{tgtSoupCntsPerGroup}.
#'
#' The different methods to determine which cells to exclude which genes for are:
#'   pCut - Exclude any cell where we can reject the null hypothesis (Poisson test) that rho<=1 for a set of genes at the \code{exCut} level.
#'   thresh - Exclude any cell for which the estimate for rho exceeds \code{exCut}.
#'
#'
#' @export 
#' @param sc A SoupChannel or SoupChannelList object.  If a SoupChannelList object, channelName must be given and must be the name of a channel.
#' @param channelName The name of a channel to use if \code{sc} is a SoupChannelList object.
#' @param nonExpressedGeneList A list containing sets of genes which can be assumed to be non-expressed in a subset of cells (see details).
#' @param useToEst A boolean matrix of dimensions length(nonExpressedGeneList) x ncol(toc) indicating which gene-sets should not be assumed to be non-expressed in each cell.  Row names must correspond to the names of \code{nonExpressedGeneList}.  If NULL, this is estimated from the data (see details).
#' @param cellGroups A vector indicating which cells to group together when estimating rho (see details).
#' @param tgtSoupCntsPerGroup When automatically constructing groups, ensure that each group has roughly this many expected soup counts.
#' @param excludeMethod Which method to use to exclude cells (see details).
#' @param exCut Cut-off used by \code{excludeMethod} (see details).
#' @return A modified version of \code{sc} with the entry \code{rhoGrouped} containing estimates of the contamination fraction for different groupings.  If \code{sc} is a SoupChannelList object the returned object has the channel named \code{channelName} modified to include the estimates of rho.
calculateContaminationFraction = function(sc,channelName,nonExpressedGeneList,useToEst=NULL,cellGroups=NULL,tgtSoupCntsPerGroup=1000,excludeMethod=c('pCut','thresh'),exCut=0.05){
  if(is(sc,'SoupChannelList')){
    if(!(channelName %in% names(sc$channels)))
      stop("sc is a SoupChannelList object, but channelName is not the name of a channel within this object")
    sc$channels[[channelName]] = calculateContaminationFraction(sc$channels[[channelName]],nonExpressedGeneList=nonExpressedGeneList,useToEst=useToEst,cellGroups=cellGroups,tgtSoupCntsPerGroup=tgtSoupCntsPerGroup,excludeMethod=excludeMethod,exCut=exCut)
    return(sc)
  }else if(!is(sc,'SoupChannel')){
    stop("sc must be a SoupChannel or SoupChannelList object")
  }
  excludeMethod = match.arg(excludeMethod)
  #Check that soup has been estimated
  if(is.null(sc$soupProfile))
    stop("Must run estimateSoup first.")
  #Convert nonExpressedGeneList to a list if we just have one set.
  if(!is.list(nonExpressedGeneList))
    nonExpressedGeneList = list(Markers=nonExpressedGeneList)
  #Collapse genes in each marker list down into one measure
  sFracs = do.call(rbind,lapply(nonExpressedGeneList,function(e) estRateLims(sum(sc$soupProfile[e,'cnts']),sum(sc$soupProfile$cnts))))
  #Do the same with cells, we don't care about other genes for estimation (other than the total counts)
  dat = do.call(rbind,lapply(nonExpressedGeneList,function(e) colSums(sc$toc[e,,drop=FALSE])))
  if(is.null(dim(dat)))
    dat = matrix(dat,nrow=1,dimnames=list(names(nonExpressedGeneList),names(dat)))
  #############################
  # Create cell mask for soup
  #For each marker set, work out which cells to exclude
  if(is.null(useToEst)){
    useToEst = t(sapply(names(nonExpressedGeneList),function(e) {
                    if(excludeMethod=='thresh'){
                      return((dat[e,]/(sc$nUMIs*sFracs[e,'est']))<exCut)
                    }else{
                      qVals = ppois(dat[e,]-1,sc$nUMIs*sFracs[e,'est'],lower.tail=FALSE)
                      qVals = p.adjust(qVals,method='BH')
                      ifelse(qVals<exCut,FALSE,TRUE)
                    }
              }))
  }
  ###############################
  # Work out how to group cells
  if(is.null(cellGroups)){
    o = order(sc$nUMIs)
    #The expected counts per cell
    expCnts = colSums(sFracs$est*useToEst)*sc$nUMIs
    #Split groups so that the expected soup count is respected
    cellGroups = cut_interval(cumsum(expCnts[o]),length=tgtSoupCntsPerGroup)
    #Reverse ordering
    cellGroups = cellGroups[match(seq_along(o),o)]
    #Drop empty levels
    cellGroups = factor(as.character(cellGroups))
    ################################
    # Deprecated bisection method.
    #Logic here is we group cells until the upper 95% confidence interval reaches required accuracy
    #tmp = colSums(sFracs$est*useToEst)
    ##Do this by recursively bisecting and finding the optimal midpoint
    #bisectGroup = function(nUMIs,tmp,lvl='',tgt){
    #  #par must be between 1 and length(nUMIs)-1
    #  penaltyFun = function(par){
    #    par = round(par)
    #    w = (seq_along(tmp)<=par)
    #    max(qbeta(0.975,1,sum(nUMIs[w]))/sum(tmp[w]),qbeta(0.975,1,sum(nUMIs[!w]))/sum(tmp[!w]))
    #  }
    #  #This will usually be faster, but I'm not 100% confident it can be used without issues so let's not use it...
    #  #fit = optim(ceiling(length(tmp)/2),penaltyFun,lower=1,upper=length(tmp)-1,method='Brent')
    #  #The exhaustive mode.  Will work, but needlessly slow for large N
    #  fit = sapply(seq(1,length(tmp)-1),penaltyFun)
    #  #Create output group mask
    #  mask = rep(lvl,length(tmp))
    #  #Test if we're done, if not split and descend a level
    #  if(min(fit)<tgt){
    #    w = seq_along(tmp)<=which.min(fit)
    #    mask[w]= bisectGroup(nUMIs[w],tmp[w],paste0(lvl,'L'),tgt=tgt)
    #    mask[!w] = bisectGroup(nUMIs[!w],tmp[!w],paste0(lvl,'R'),tgt=tgt)
    #  }
    #  return(mask)
    #}
    ##Group things at required power threshold within bins of similar nUMIs
    #grps = bisectGroup(nUMIs[o],tmp[o],'',tgt=tgtAccuracy)
    #Re-order them in the input way
    #cellGroups = grps[match(seq_along(o),o)]
    ##################################
  }
  ##############
  # Estimation
  obsSoupCnts = estRateLims(colSums(dat*useToEst),sc$nUMIs)*sc$nUMIs
  expSoupCnts = colSums(sFracs$est*useToEst)*sc$nUMIs
  #Global estimate
  globRho = estGrouped(obsSoupCnts$est,expSoupCnts,sc$nUMIs,rep('Global',length(sc$nUMIs)))
  if('Global' %in% cellGroups)
    stop('Cell grouping named "Global" is reserved.  Please rename cell group.')
  # calculate rho in groups
  groupRhos = estGrouped(obsSoupCnts$est,expSoupCnts,sc$nUMIs,cellGroups)
  sc$rhoGrouped = (rbind(globRho,groupRhos))
  return(sc)
}
